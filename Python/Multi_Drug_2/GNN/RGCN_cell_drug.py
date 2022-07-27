import argparse
import torch
import torch.nn.functional as F
from torch.nn import Linear

import torch_geometric.transforms as T

from Create_Het_Graph import MyGraphDataset
from torch_geometric.nn import RGCNConv
#from sklearn.metrics import roc_auc_score
from torch_geometric.nn import SAGEConv, to_hetero


# load the dataset
dataset = MyGraphDataset(root= "Raw_Data_From_R/Cellline_drug_Net")
print(dataset.data)
print(dataset.data.edge_types)
print(dataset.data.node_types)

data = dataset.data



parser = argparse.ArgumentParser()
parser.add_argument('--use_weighted_loss', action='store_true',
                    help='Whether to use weighted MSE loss.')
args = parser.parse_args()



# Add user node features for message passing:
#data['user'].x = torch.eye(data['user'].num_nodes)
#del data['cellline'].num_nodes

# Add a reverse ('movie', 'rev_rates', 'user') relation for message passing:
data = T.ToUndirected()(data)

# Perform a link-level split into training, validation, and test edges:
train_data, val_data, test_data = T.RandomLinkSplit(
    num_val=0.1,
    num_test=0.1,
    neg_sampling_ratio=1,
    add_negative_train_samples = True,
    edge_types=[('cellline', 'sen', 'drug')],
    rev_edge_types=[('drug', 'rev_sen', 'cellline')]
)(data)


def dict_conversion(data):
    s = data.metadata()[1]
    edge_index = torch.cat([data[s[0]].edge_index, 
                         data[s[1]].edge_index,
                         data[s[2]].edge_index, 
                         data[s[3]].edge_index], dim=-1)
    
    E = edge_index.shape[1]
    
    edge_type = torch.zeros(E, dtype=torch.float)
    start = 0
    for i in range(4):
        end = start + data[s[i]].edge_index.size(1) 
        edge_type[start : end] = i
        start = end + 1
    

    
    
# We have an unbalanced dataset with many labels for rating 3 and 4, and very

class GNNEncoder(torch.nn.Module):
    def __init__(self, hidden_channels, out_channels = 32, num_relations = 4):
        super().__init__()
        self.conv1 = RGCNConv((3831,1024), hidden_channels, num_relations = num_relations)
        self.conv2 = RGCNConv((3831,1024), out_channels, num_relations = num_relations)

    def forward(self, x, edge_index, edge_type = [('cellline','sen','drug')]):
        x = self.conv1(x, edge_index, edge_type = edge_type).relu()
        x = self.conv2(x, edge_index, edge_type = edge_type)
        return x


class EdgeDecoder(torch.nn.Module):
    def __init__(self, hidden_channels):
        super().__init__()
        self.lin1 = Linear(2 * hidden_channels, hidden_channels)
        self.lin2 = Linear(hidden_channels, 1)

    def forward(self, z_dict, edge_label_index):
        row, col = edge_label_index
        z = torch.cat([z_dict['cellline'][row], z_dict['drug'][col]], dim=-1)

        z = self.lin1(z).relu()
        z = self.lin2(z)
        return z.view(-1)
    

class Model(torch.nn.Module):
    def __init__(self, hidden_channels):
        super().__init__()
        self.encoder = GNNEncoder(hidden_channels, hidden_channels,num_relations = 3)
        self.encoder = to_hetero(self.encoder, data.metadata(), aggr='sum')
        self.decoder = EdgeDecoder(hidden_channels)

    def forward(self, x_dict, edge_index_dict, edge_label_index):
        z_dict = self.encoder(x_dict, edge_index_dict)
        return self.decoder(z_dict, edge_label_index)

# main:

model = Model(hidden_channels=32)

# Due to lazy initialization, we need to run one model step so the number
# of parameters can be inferred:
with torch.no_grad():
    model.encoder(train_data.x_dict, train_data.edge_index_dict)

optimizer = torch.optim.Adam(model.parameters(), lr=0.01)


def train():
    model.train()
       
    optimizer.zero_grad()
    pred = model(train_data.x_dict, train_data.edge_index_dict,
                 train_data['cellline', 'sen', 'drug'].edge_label_index)
    target = train_data['cellline', 'sen', 'drug'].edge_label
    loss = F.binary_cross_entropy_with_logits(pred, target)
    loss.backward()
    optimizer.step()
    return float(loss)


@torch.no_grad()
def test(d):
    model.eval()
    pred = model(d.x_dict, d.edge_index_dict,
                 d['cellline', 'sen', 'drug'].edge_label_index)
    target = d['cellline', 'sen', 'drug'].edge_label.float()
    loss = F.binary_cross_entropy_with_logits(pred, target)
    return float(loss)


for epoch in range(1, 10):
    loss = train()
    train_rmse = test(train_data)
    val_rmse = test(val_data)
    test_rmse = test(test_data)
    print(f'Epoch: {epoch:03d}, Loss: {loss:.4f}, Train: {train_rmse:.4f}, '
          f'Val: {val_rmse:.4f}, Test: {test_rmse:.4f}')

