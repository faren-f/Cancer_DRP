#branch:Loss_NA
#import os
#os.chdir("Desktop/Cancer_DRP/Python/GNN")

import torch
import numpy as np
import torch_geometric.transforms as T

from torch_geometric.nn import SAGEConv, to_hetero, GraphConv
#from numpy import savetxt
from RandomLinkSplit_modified_NewNode import RandomLinkSplit
from Create_Het_Graph_Regression import MyGraphDataset

# load the dataset
dataset = MyGraphDataset(root= "Raw_Data_From_R/Cellline_drug_Net/Regression")
print(dataset.data)
print(dataset.data.edge_types)
print(dataset.data.node_types)



data = dataset.data

# Add a reverse ('item', 'rev_rates', 'user') relation for message passing:
data = T.ToUndirected()(data)
del data['drug', 'rev_sen', 'cellline'].edge_label  # Remove "reverse" label.

# Perform a link-level split into training, validation, and test edges:
splitdata = RandomLinkSplit(
    num_val=0.1,
    num_test=0.1,
    is_undirected = False,
    neg_sampling_ratio=0,
    edge_types=[('cellline', 'sen', 'drug')],
    rev_edge_types=[('drug', 'rev_sen', 'cellline')])
train_data, val_data, test_data = splitdata(data)

print(train_data)
print(val_data)
print(test_data)



def weighted_mse_loss(pred, target, weight=None):
    weight = 1. if weight is None else weight[target].to(pred.dtype)
    return (weight * (pred - target.to(pred.dtype)).pow(2)).mean()


class GNNEncoder(torch.nn.Module):
    def __init__(self, hidden_channels, out_channels):
        super().__init__()
        self.conv1 = GraphConv((-1, -1), hidden_channels)
        #self.conv2 = SAGEConv((-1, -1), hidden_channels)
        self.conv2 = GraphConv((-1, -1), out_channels)

    def forward(self, x, edge_index, edge_weight):
        x = self.conv1(x, edge_index,  edge_weight).relu()
        #x = self.conv2(x, edge_index).relu()
        x = self.conv2(x, edge_index,  edge_weight)

        return x
    

class Model(torch.nn.Module):
    def __init__(self, hidden_channels):
        super().__init__()
        self.encoder = GNNEncoder(hidden_channels, hidden_channels)
        self.encoder = to_hetero(self.encoder, data.metadata(), aggr='sum')
        self.sigmoid = torch.nn.Sigmoid()

    def forward(self, x_dict, edge_index_dict, edge_weight_dict, edge_label_index):
        z_dict = self.encoder(x_dict, edge_index_dict, edge_weight_dict)
        row1, row2 = edge_label_index
        inner_product = (z_dict['cellline'][row1] * z_dict['drug'][row2]).sum(dim=-1)  # dot product 
        #sigmoid_inner_product = self.sigmoid(inner_product)
        return inner_product
    # This function gives the embedding space that we need for interpretation
    def forward_encoder(self, x_dict, edge_index_dict, edge_weight_dict):
        z_dict = self.encoder(x_dict, edge_index_dict, edge_weight_dict)
        return z_dict


# main:

model = Model(hidden_channels=80)


# Due to lazy initialization, we need to run one model step so the number
# of parameters can be inferred:
with torch.no_grad():
    model.encoder(train_data.x_dict, train_data.edge_index_dict, train_data.edge_weight_dict)

optimizer = torch.optim.Adam(model.parameters(), lr=0.0006)


def train():
    model.train()
    optimizer.zero_grad()
    pred = model(train_data.x_dict, train_data.edge_index_dict,
                 train_data.edge_weight_dict,
                 train_data['cellline', 'sen', 'drug'].edge_label_index)
    target = train_data['cellline', 'sen', 'drug'].edge_label
    
    loss = weighted_mse_loss(pred, target, weight= None)
    loss.backward()
    optimizer.step()
    return float(loss)


@torch.no_grad()
def test(d):
    model.eval()
    pred = model(d.x_dict, d.edge_index_dict,d.edge_weight_dict,
                 d['cellline', 'sen', 'drug'].edge_label_index)
    target = d['cellline', 'sen', 'drug'].edge_label.float()
    # rmse = F.mse_loss(pred, target).sqrt()
    cor = np.corrcoef(pred, target)[0,1]
    return float(cor)    

@torch.no_grad()
def embbeding(d):
    model.eval()
    z = model.forward_encoder(d.x_dict, d.edge_index_dict, d.edge_weight_dict)
    return (z)  
 

for epoch in range(1, 60):
    loss = train()
    train_cor = test(train_data)
    val_cor = test(val_data) 
    test_cor = test(test_data)
    print(f'Epoch: {epoch:03d}, Loss: {loss:.4f}, Train: {train_cor:.4f}, '
          f'Val: {val_cor:.4f}, Test: {test_cor:.4f}')        


z = embbeding(data)
Embedding_cellline = z["cellline"].detach().numpy()
Embedding_drug = z["drug"].detach().numpy()

#savetxt('Saved_data/RGCN/Embedding_cellline_prod.csv', Embedding_cellline,delimiter=',')
#savetxt('Saved_data/RGCN/Embedding_drug_prod_.csv', Embedding_drug,delimiter=',') 
