from Create_Graph import MyGraphDataset
import torch
import torch.nn.functional as F
from torch_geometric.nn import GCNConv
from sklearn import metrics


class GCN(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.conv1 = GCNConv(dataset.num_node_features, 20)
        #self.conv2 = GCNConv(10, 5)
        self.conv2 = GCNConv(20, 2)

    def forward(self, data):
        x, edge_index = data.x, data.edge_index

        x = self.conv1(x, edge_index)
        x = F.leaky_relu(x)
        x = F.dropout(x, training=self.training)
        x = self.conv2(x, edge_index)
        #x = F.leaky_relu(x)
        #x = F.dropout(x, training=self.training)
        #x = self.conv3(x, edge_index)


        return F.log_softmax(x, dim=1)
    
    
## main

dataset = MyGraphDataset(root= "Raw_Data_From_R/")
model = GCN()
data = dataset.data


optimizer = torch.optim.Adam(model.parameters(), lr=0.01, weight_decay=5e-4)

model.train()
for epoch in range(100):
    print(epoch)
    optimizer.zero_grad()
    out = model(data)
    loss = F.nll_loss(out, data.y)
    #loss = F.mse_loss(out, data.y)
    loss.backward()
    optimizer.step()
    
## Evaluation
model.eval()
pred = model(data).argmax(dim=1)
#correct = (pred == data.y).sum()

auc = metrics.roc_auc_score(data.y, pred)
print(auc)


#acc = int(correct) / dataset.N_Nodes
#print(f'Accuracy: {acc:.4f}')



