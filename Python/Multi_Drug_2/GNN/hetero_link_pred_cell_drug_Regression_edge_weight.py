import argparse
import torch
import torch.nn.functional as F
from torch.nn import Linear

import torch_geometric.transforms as T
#from torch_geometric.transforms import RandomLinkSplit
from RandomLinkSplit_modified_NewNode import RandomLinkSplit

from Create_Het_Graph_Regression import MyGraphDataset
from torch_geometric.nn import SAGEConv, to_hetero, LEConv, GraphConv
import numpy as np
#from numpy import savetxt

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
splitdata = RandomLinkSplit(
    num_val=0.1,
    num_test=0.1,
    neg_sampling_ratio=0,
    edge_types=[('cellline', 'sen', 'drug')],
    rev_edge_types=[('drug', 'rev_sen', 'cellline')])
    
train_data, val_data, test_data = splitdata(data) 


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
        self.conv1 = LEConv((-1, -1), hidden_channels)
        #self.conv2 = LEConv((-1, -1), hidden_channels)
        #self.conv2 = LEConv((-1, -1), out_channels)

    def forward(self, x, edge_index, edge_weight):
        x = self.conv1(x, edge_index, edge_weight)
        #x = self.conv2(x, edge_index, edge_weight)
        #x = self.conv3(x, edge_index, edge_weight)

        return x


class EdgeDecoder(torch.nn.Module):
    def __init__(self, hidden_channels):
        super().__init__()
        self.lin1 = Linear(2 * hidden_channels, hidden_channels)
        #self.lin2 = Linear(hidden_channels, hidden_channels)
        self.lin2 = Linear(hidden_channels, 1)

    def forward(self, z_dict, edge_label_index):
        row1, row2 = edge_label_index
        z = torch.cat([z_dict['cellline'][row1], z_dict['drug'][row2]], dim=-1)

        z = self.lin1(z).sigmoid()
        z = self.lin2(z).sigmoid()
        #z = self.lin3(z).sigmoid()
        return z.view(-1)
    

class Model(torch.nn.Module):
    def __init__(self, hidden_channels):
        super().__init__()
        self.encoder = GNNEncoder(hidden_channels, hidden_channels)
        self.encoder = to_hetero(self.encoder, data.metadata(), aggr='mean')
        self.decoder = EdgeDecoder(hidden_channels)

    def forward(self, x_dict, edge_index_dict, edge_weight_dict, edge_label_index):
        z_dict = self.encoder(x_dict, edge_index_dict, edge_weight_dict)
        return self.decoder(z_dict, edge_label_index)
    
    def forward_encoder(self, x_dict, edge_index_dict, edge_weight_dict):
        z_dict = self.encoder(x_dict, edge_index_dict, edge_weight_dict)
        return z_dict
# main:

Test_nodes = splitdata.test_nodes.numpy()
Val_nodes = splitdata.val_nodes.numpy()
Train_nodes = splitdata.train_nodes.numpy()

sen = dataset.sen
sen_Test = sen[Test_nodes,] 
sen_Val = sen[Val_nodes,] 
sen_Train = sen[Train_nodes,] 


model = Model(hidden_channels=100)


# Due to lazy initialization, we need to run one model step so the number
# of parameters can be inferred:
with torch.no_grad():
    model.encoder(train_data.x_dict, train_data.edge_index_dict, train_data.edge_weight_dict)

optimizer = torch.optim.Adam(model.parameters(), lr=0.0007)


def train():
    model.train()
    optimizer.zero_grad()
    pred = model(train_data.x_dict, train_data.edge_index_dict,
                 train_data.edge_weight_dict,
                 train_data['cellline', 'sen', 'drug'].edge_label_index)
    target = train_data['cellline', 'sen', 'drug'].edge_label
    
    loss = weighted_mse_loss(pred, target, weight)
    loss.backward()
    optimizer.step()
    return float(loss)


@torch.no_grad()
def test(data,sen,sen_data):
    model.eval()
    pred = model(data.x_dict, data.edge_index_dict, data.edge_weight_dict,
                 data['cellline', 'sen', 'drug'].edge_label_index)
    
    
    prediction = np.zeros([sen_data.shape[0],sen_data.shape[1]])
    
    for t in range(0,sen_data.shape[0]):
        
        N = sen.shape[1] - sum(np.isnan(sen_data[t,:]))
        Ind_nan = np.where(np.isnan(sen_data[t,:]))
        Pred = pred[0:N].detach().numpy()
        prediction[t,Ind_nan] = np.nan
        prediction[t,~np.isnan(prediction[t,])] = Pred
        pred = pred[N::]
        
    cor = np.array([ ])
    loss = np.array([ ])    
    for j in range(0, sen_data.shape[1]):
        Target = sen_data[~np.isnan(sen_data[:,j]),j]
        pred = prediction[~np.isnan(prediction[:,j]),j]
        #Target = sen_data[j, ~np.isnan(sen_data[j,:])]
        #pred = prediction[j, ~np.isnan(prediction[j,:])]
        cor = np.append(cor, np.corrcoef(pred, Target)[0,1])
        loss = np.append(loss, F.mse_loss(torch.tensor(pred), torch.tensor(Target)).sqrt())

    return (np.nanmean(cor),float(np.nanmean(loss)))
        
              
#loss_test = np.append(loss_test, weighted_mse_loss(Pred, Target, weight))
# print(cor[t])

#target = d['cellline', 'sen', 'drug'].edge_label.float()
#rmse = F.mse_loss(pred, target).sqrt()
#cor = np.corrcoef(pred, target)[0,1]
#loss_test = weighted_mse_loss(pred, target, weight)



@torch.no_grad()
def embbeding(d):
    model.eval()
    z = model.forward_encoder(d.x_dict, d.edge_index_dict, d.edge_weight_dict)
    return (z)  

    

for epoch in range(1, 100):
    print(f"Epoch: {epoch+1}")
    loss = train()
    train_corr,train_loss = test(train_data, sen, sen_Train)
    val_corr,val_loss = test(val_data, sen, sen_Val) 
    test_corr,test_loss = test(test_data, sen, sen_Test)
    print(f'Epoch: {epoch:03d}, Loss: {loss:.4f}, Train: {train_corr:.4f}, '
          f'Val: {val_corr:.4f}, Test: {test_corr:.4f}') 
       
    
        

#z = embbeding(data)
#Embedding_cellline = z["cellline"].detach().numpy()
#Embedding_drug = z["drug"].detach().numpy()

#savetxt('Saved_data/RGCN/Embedding_cellline.csv', Embedding_cellline,delimiter=',')
#savetxt('Saved_data/RGCN/Embedding_drug.csv', Embedding_drug,delimiter=',') 




