import argparse
import torch
import torch.nn.functional as F
from torch.nn import Linear

import torch_geometric.transforms as T
#from torch_geometric.transforms import RandomLinkSplit
from RandomLinkSplit_modified_NewNode import RandomLinkSplit

from Create_Het_Graph_Regression import MyGraphDataset
from torch_geometric.nn import SAGEConv, to_hetero
import numpy as np
from numpy import savetxt

# load the dataset
dataset = MyGraphDataset(root= "Raw_Data_From_R/Cellline_drug_Net/Regression")
print(dataset.data)
print(dataset.data.edge_types)
print(dataset.data.node_types)

data = dataset.data



parser = argparse.ArgumentParser()
parser.add_argument('--use_weighted_loss', action='store_true',
                    help='Whether to use weighted MSE loss.')
args = parser.parse_args()


# Add a reverse ('movie', 'rev_rates', 'user') relation for message passing:
data = T.ToUndirected()(data)
del data['drug', 'rev_sen', 'cellline'].edge_label  # Remove "reverse" label.

# Perform a link-level split into training, validation, and test edges:
train_data, val_data, test_data = RandomLinkSplit(
    num_val=0.1,
    num_test=0.1,
    neg_sampling_ratio=0,
    edge_types=[('cellline', 'sen', 'drug')],
    rev_edge_types=[('drug', 'rev_sen', 'cellline')],
)(data)


# We have an unbalanced dataset with many labels for rating 3 and 4, and very
# few for 0 and 1. Therefore we use a weighted MSE loss.
if args.use_weighted_loss:
    weight = torch.bincount(train_data['cellline','sen', 'drug'].edge_label)
    weight = weight.max() / weight
else:
    weight = None


def weighted_mse_loss(pred, target, weight=None):
    weight = 1. if weight is None else weight[target].to(pred.dtype)
    return (weight * (pred - target.to(pred.dtype)).pow(2)).mean()


class GNNEncoder(torch.nn.Module):
    def __init__(self, hidden_channels, out_channels):
        super().__init__()
        self.conv1 = SAGEConv((-1, -1), hidden_channels)
        #self.conv2 = SAGEConv((-1, -1), hidden_channels)
        self.conv2 = SAGEConv((-1, -1), out_channels)

    def forward(self, x, edge_index):
        x = self.conv1(x, edge_index).relu()
        #x = self.conv2(x, edge_index).relu()
        x = self.conv2(x, edge_index)

        return x


class EdgeDecoder(torch.nn.Module):
    def __init__(self, hidden_channels):
        super().__init__()
        self.lin1 = Linear(2 * hidden_channels, hidden_channels)
        self.lin2 = Linear(hidden_channels, 1)

    def forward(self, z_dict, edge_label_index):
        row1, row2 = edge_label_index
        z = torch.cat([z_dict['cellline'][row1], z_dict['drug'][row2]], dim=-1)

        z = self.lin1(z).relu()
        z = self.lin2(z)
        return z.view(-1)
    

class Model(torch.nn.Module):
    def __init__(self, hidden_channels):
        super().__init__()
        self.encoder = GNNEncoder(hidden_channels, hidden_channels)
        self.encoder = to_hetero(self.encoder, data.metadata(), aggr='sum')
        self.decoder = EdgeDecoder(hidden_channels)

    def forward(self, x_dict, edge_index_dict, edge_label_index):
        z_dict = self.encoder(x_dict, edge_index_dict)
        return self.decoder(z_dict, edge_label_index)
    
    def forward_encoder(self, x_dict, edge_index_dict):
        z_dict = self.encoder(x_dict, edge_index_dict)
        return z_dict
# main:

model = Model(hidden_channels=100)


# Due to lazy initialization, we need to run one model step so the number
# of parameters can be inferred:
with torch.no_grad():
    model.encoder(train_data.x_dict, train_data.edge_index_dict)

optimizer = torch.optim.Adam(model.parameters(), lr=0.0007)


def train():
    model.train()
    optimizer.zero_grad()
    pred = model(train_data.x_dict, train_data.edge_index_dict,
                 train_data['cellline', 'sen', 'drug'].edge_label_index)
    target = train_data['cellline', 'sen', 'drug'].edge_label
    
    loss = weighted_mse_loss(pred, target, weight)
    loss.backward()
    optimizer.step()
    return float(loss)


@torch.no_grad()
def test(d):
    model.eval()
    pred = model(d.x_dict, d.edge_index_dict,
                 d['cellline', 'sen', 'drug'].edge_label_index)
    target = d['cellline', 'sen', 'drug'].edge_label.float()
    # rmse = F.mse_loss(pred, target).sqrt()
    cor = np.corrcoef(pred, target)[0,1]
    loss_test = weighted_mse_loss(pred, target, weight)

    return float(cor)    

@torch.no_grad()
def embbeding(d):
    model.eval()
    z = model.forward_encoder(d.x_dict, d.edge_index_dict)
    return (z)  
 

for epoch in range(1, 60):
    loss = train()
    train_corr = test(train_data)
    val_corr = test(val_data) 
    test_corr = test(test_data)
    print(f'Epoch: {epoch:03d}, Loss: {loss:.4f}, Train: {train_corr:.4f}, '
          f'Val: {val_corr:.4f}, Test: {test_corr:.4f}')        


#z = embbeding(data)
#Embedding_cellline = z["cellline"].detach().numpy()
#Embedding_drug = z["drug"].detach().numpy()

#savetxt('Saved_data/RGCN/Embedding_cellline.csv', Embedding_cellline,delimiter=',')
#savetxt('Saved_data/RGCN/Embedding_drug.csv', Embedding_drug,delimiter=',') 
