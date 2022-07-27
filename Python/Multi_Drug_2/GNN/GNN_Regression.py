from Create_Graph import MyGraphDataset
import torch.nn as nn
import torch.nn.functional as F
import torch_geometric.nn as pyg_nn
import torch.optim as optim
import numpy as np
import torch
import matplotlib.pyplot as plt

## class definition

class GNNStack(nn.Module):
    def __init__(self, hidden_dim_GNN, hidden_dim_MLP, output_dim, num_layers):
        super(GNNStack, self).__init__()
        self.num_layers = num_layers
        input_dim = hidden_dim_GNN[0]
        self.convs = nn.ModuleList()
        self.convs.append(self.build_conv_model(input_dim, hidden_dim_GNN[1]))
        self.lns = nn.ModuleList()

        for l in range(num_layers-1):
            self.convs.append(self.build_conv_model(hidden_dim_GNN[l+1], hidden_dim_GNN[l+2]))
            self.lns.append(nn.LayerNorm(hidden_dim_GNN[l+1]))

        # post-message-passing
        self.dropout = 0.2
        
        self.post_mp = nn.Sequential(
            nn.Linear(hidden_dim_GNN[-1], hidden_dim_MLP[0]), nn.Dropout(self.dropout), nn.ReLU(),  
            #nn.Linear(hidden_dim_MLP[0], hidden_dim_MLP[1]), nn.Dropout(self.dropout), nn.ReLU(),
            nn.Linear(hidden_dim_MLP[0], output_dim))


    def build_conv_model(self, input_dim, hidden_dim):
        # refer to pytorch geometric nn module for different implementation of GNNs.
        #return pyg_nn.GINConv(nn.Sequential(nn.Linear(input_dim, hidden_dim),
        #                          nn.LeakyReLU(), nn.Linear(hidden_dim, hidden_dim)))
        #return pyg_nn.GATConv(input_dim, hidden_dim)
        #return pyg_nn.GCNConv(input_dim, hidden_dim)
        return pyg_nn.SAGEConv(input_dim, hidden_dim)
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

        pred = self.post_mp(x)

        return emb, pred

    def loss(self, pred, label):
        return F.mse_loss(pred, label)
    

## main

dataset = MyGraphDataset(root= "Raw_Data_From_R/")
data = dataset.data

    
#build model
model = GNNStack([dataset.num_node_features,50,5,1],[50], 1, num_layers = 3)
opt = optim.Adam(model.parameters(), lr=0.001)

## Balance the number of nodes
label = data.y
Cor = np.array([])

for rep in range(1):
    
    # Seperate train and test sets
    num_train_nodes = round(0.8 * dataset.N_Nodes)
    perm = torch.randperm(dataset.N_Nodes)
    train_idx = perm[:num_train_nodes]
    test_idx = perm[num_train_nodes:]
    Loss = np.array([])      
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
        Loss = np.append(Loss,loss.detach().numpy())

         # gradients = backward pass
        loss.backward()
         # opdate weights
        opt.step()
    
    with torch.no_grad():
        emb, pred = model(data)        

    cor = np.corrcoef(data.y[test_idx], np.squeeze(pred[test_idx]))
    print(cor[0,1])

    Cor = np.append(Cor,cor[0,1])


print("Final Result: ")
print("Corelation mean = ", Cor.mean())
print("Corelation sd = ", Cor.std())

plt.plot(Loss)







    
       