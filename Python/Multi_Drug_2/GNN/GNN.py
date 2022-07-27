from Create_Graph import MyGraphDataset
import torch.nn as nn
import torch.nn.functional as F
import torch_geometric.nn as pyg_nn
import torch.optim as optim
from sklearn import metrics
import numpy as np
from itertools import chain
from random import sample
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
        self.dropout = 0
        self.post_mp = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim), nn.Dropout(self.dropout), nn.ReLU(),   
            # nn.Linear(hidden_dim, hidden_dim), nn.Dropout(self.dropout), nn.ReLU(),
            nn.Linear(hidden_dim, output_dim))


    def build_conv_model(self, input_dim, hidden_dim):
        # refer to pytorch geometric nn module for different implementation of GNNs.
        #return pyg_nn.GINConv(nn.Sequential(nn.Linear(input_dim, hidden_dim),
        #                          nn.LeakyReLU(), nn.Linear(hidden_dim, hidden_dim)))
        #return pyg_nn.GATConv(input_dim, hidden_dim)
        return pyg_nn.GCNConv(input_dim, hidden_dim)
        #return pyg_nn.SAGEConv(input_dim, hidden_dim)
        #return pyg_nn.GraphConv(input_dim, hidden_dim)
        #return pyg_nn.LEConv(input_dim, hidden_dim)
        #return pyg_nn.GCN2Conv(channels = 1, alpha=1)




    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        for i in range(self.num_layers):
            x = self.convs[i](x, edge_index = edge_index)
            emb = x
            x = F.relu(x)
            x = F.dropout(x, p=self.dropout, training=self.training)
            #if not i == self.num_layers - 1:
               #x = self.lns[i](x)

        x = self.post_mp(x)

        return emb, F.log_softmax(x, dim=1)

    def loss(self, pred, label):
        return F.nll_loss(pred, label)
    

## main

dataset = MyGraphDataset(root= "Raw_Data_From_R/")
data = dataset.data

    
#build model
model = GNNStack(dataset.num_node_features, 100, 2, num_layers = 4)
opt = optim.Adam(model.parameters(), lr=0.001)

## Balance the number of nodes
label = data.y
AUC = np.array([])

for rep in range(1):
    
    # Under sample the majarity class
    n1 = sum(label==0)
    n2 = sum(label==1)
    # Number of data that should be remove (for under-sampling)
    m = abs(n2-n1)
    
    # gives the indeces of majarity class 
    if n1>n2:
        ind_majarity = np.where(label==0)
    else:   
        ind_majarity = np.where(label==1)
        
    
    ## Turn tuple of tuples to list because "sample" accepts sequence in its first argument 
    ind_majarity = list(chain(*ind_majarity)) 
    # indeces that we want to remove from the majarity class
    ind_remove = sample((ind_majarity),k=m)
    
    all_ind = range(len(data.y))
    ind_to_keep = np.delete(all_ind, ind_remove)
    
    # Seperate train and test sets
    num_train_nodes = round(0.8 * len(ind_to_keep))
    perm = torch.randperm(len(ind_to_keep))
    train_idx = perm[:num_train_nodes]
    test_idx = perm[num_train_nodes:]
          
    # train
    for epoch in range(50):
        print(epoch)
        total_loss = 0
        model.train()
                    
        # zero gradients
        opt.zero_grad()
        # pridiction = forward pass
        emb, pred = model(data)
                    
         
         # Loss
        loss = model.loss(pred[train_idx], label[train_idx])
        #loss = model.loss(pred, label)
         # gradients = backward pass
        loss.backward()
         # opdate weights
        opt.step()
    
    with torch.no_grad():
        emb, pred = model(data)        
        pred = pred.argmax(dim=1)

    auc = metrics.roc_auc_score(data.y[test_idx], pred[test_idx])
    print(auc)

    AUC = np.append(AUC,auc)


print("Final Result: ")
print("AUC mean = ", AUC.mean())
print("AUC sd = ", AUC.std())










    
       