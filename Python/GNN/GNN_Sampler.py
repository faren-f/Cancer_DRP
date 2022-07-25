from Create_Graph import MyGraphDataset
import torch.nn as nn
import torch.nn.functional as F
import torch_geometric.nn as pyg_nn
import torch.optim as optim
from sklearn import metrics
from torch.utils.data import WeightedRandomSampler
import numpy as np
import torch.utils.data
from torch_geometric.data import DataLoader
import torch


## class definition

class GNNStack(nn.Module):
    def __init__(self, input_dim, hidden_dim, output_dim, num_layers=3):
        super(GNNStack, self).__init__()
        self.num_layers = num_layers
        self.convs = nn.ModuleList()
        self.convs.append(self.build_conv_model(input_dim, hidden_dim))
        self.lns = nn.ModuleList()

        for l in range(num_layers-1):
            self.convs.append(self.build_conv_model(hidden_dim, hidden_dim))
            self.lns.append(nn.LayerNorm(hidden_dim))

        # post-message-passing
        self.dropout = .25
        self.post_mp = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim), nn.Dropout(self.dropout), nn.ReLU(),   
            #nn.Linear(hidden_dim, hidden_dim), nn.Dropout(self.dropout),
            nn.Linear(hidden_dim, output_dim))


    def build_conv_model(self, input_dim, hidden_dim):
        # refer to pytorch geometric nn module for different implementation of GNNs.
        #return pyg_nn.GINConv(nn.Sequential(nn.Linear(input_dim, hidden_dim),
        #                          nn.LeakyReLU(), nn.Linear(hidden_dim, hidden_dim)))
        #return pyg_nn.GATConv(input_dim, hidden_dim)
        #return pyg_nn.GCNConv(input_dim, hidden_dim)
        #return pyg_nn.SAGEConv(input_dim, hidden_dim)
        return pyg_nn.GraphConv(input_dim, hidden_dim)
        #return pyg_nn.LEConv(input_dim, hidden_dim)
        #return pyg_nn.GCN2Conv(channels = 1, alpha=1)




    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        for i in range(self.num_layers):
            x = self.convs[i](x, edge_index = edge_index)
            emb = x
            #x = F.relu(x)
            #x = F.dropout(x, p=self.dropout, training=self.training)
            #if not i == self.num_layers - 1:
               #x = self.lns[i](x)

        x = self.post_mp(x)

        return emb, F.log_softmax(x, dim=1)

    def loss(self, pred, label):
        return F.nll_loss(pred, label)
    




## main

dataset = MyGraphDataset(root= "Raw_Data_From_R/")
data = dataset.data
data_size = dataset.num_node_features
class_sample_count = np.unique(dataset.node_label, return_counts=True)[1]
weight = 1. / class_sample_count
samples_weight = np.array([weight[t] for t in dataset.node_label])
samples_weight = torch.from_numpy(samples_weight)
sampler = WeightedRandomSampler(samples_weight, len(samples_weight))

loader = DataLoader(dataset.data, batch_size=10)


    
#build model
model = GNNStack(dataset.num_node_features, 1, 2, num_layers = 3)
opt = optim.Adam(model.parameters(), lr=0.01)
      
# train
for epoch in range(50):
    print(epoch)
    total_loss = 0
    model.train()
    for batch in loader:

                
        # zero gradients
        opt.zero_grad()
        # pridiction = forward pass
        embedding, pred = model(batch)
                
        # seprate labels from batch that is part of x
        label = batch.y
        # Loss
        loss = model.loss(pred, label)
        # gradients = backward pass
        loss.backward()
        # opdate weights
        opt.step()


embedding, pred = model(data)

## Evaluation

pred = pred.argmax(dim=1)
#correct = (pred == data.y).sum()

auc = metrics.roc_auc_score(data.y, pred)
print(auc)








    
       




            
            
            
            
            
            
                
            
            
            
            
            
            
           