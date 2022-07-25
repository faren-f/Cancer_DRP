import torch
import torch.nn.functional as F
from torch.nn import Linear

import torch_geometric.transforms as T
#from torch_geometric.transforms import RandomLinkSplit

from Create_Het_Graph_Classifier import MyGraphDataset
from torch_geometric.nn import SAGEConv, to_hetero
from torch_geometric.nn import GCNConv
from torch_geometric.utils import negative_sampling
from sklearn import metrics
import numpy as np


# load the dataset
dataset = MyGraphDataset(root= "Raw_Data_From_R/Cellline_drug_Net/Classifier")
print(dataset.data)
print(dataset.data.edge_types)
print(dataset.data.node_types)

data = dataset.data


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
    rev_edge_types=[('drug', 'rev_sen', 'cellline')],
)(data)


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
        self.lin1 = Linear(2 * hidden_channels, 1)
        #self.lin2 = Linear(hidden_channels, 1)

    def forward(self, z_dict, edge_label_index):
        row, col = edge_label_index
        z = torch.cat([z_dict['cellline'][row], z_dict['drug'][col]], dim=-1)

        # z = self.lin1(z).relu()
        z = self.lin1(z)
        #z = self.lin2(z)
        return z.view(-1)      # view(-1) makes it flattened
    

class Model(torch.nn.Module):
    def __init__(self, hidden_channels):
        super().__init__()
        self.encoder = GNNEncoder(hidden_channels, hidden_channels)
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

optimizer = torch.optim.Adam(model.parameters(), lr=0.001)


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
    for i in range(len(pred)):
        if pred[i]>=0:
            pred[i]=1
        else:
            pred[i]=0
        
        
        
    target = d['cellline', 'sen', 'drug'].edge_label.float()
    
    loss = F.binary_cross_entropy_with_logits(pred, target)
    return float(loss)


for epoch in range(1, 30):
    loss = train()
    #train_loss = test(train_data)
    #val_loss = test(val_data)
    test_loss = test(test_data)
    #print(f'Epoch: {epoch:03d}, Loss: {loss:.4f}, Train: {train_loss:.4f}, '
          #f'Val: {val_loss:.4f}, Test: {test_loss:.4f}')
    print(f'Epoch: {epoch:03d}, Loss: {loss:.4f}, Test: {test_loss:.4f}')      
    print(f'Epoch: {epoch:03d}, Loss: {loss:.4f}')
                
      


# Final calculation of AUC
#def testfinal(d):
    #model.eval()
    #pred = model(d.x_dict, d.edge_index_dict,
                 #d['cellline', 'sen', 'drug'].edge_label_index)
    #for i in range(len(pred)):
        #if pred[i]>0:
            #pred[i]=1
        #else:
            #pred[i]=0
        
    #target = d['cellline', 'sen', 'drug'].edge_label.float()
    
    
with torch.no_grad():
    auc = metrics.roc_auc_score(pred, target)
    acc = metrics.accuracy_score(pred,target)
    sen = metrics.recall_score(pred,target)
    c  = np.corrcoef(pred, target)
    confusion = metrics.confusion_matrix(pred, target)
print(c)
print(auc)
print(acc)
print(sen)
print(confusion)


#train_auc = testfinal(train_data)
#val_auc = testfinal(val_data)
#test_auc = testfinal(test_data)
#print(f'Train: {train_auc:.4f}, 'f'Val: {val_auc:.4f}, Test: {test_auc:.4f}')