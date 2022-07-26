import numpy as np

import torch
from torch import nn

from torch.utils.data import DataLoader
from Create_CustomDataset import CustomDataset


class MLP(nn.Module):
    def __init__(self, input_size, hidden_size, output_size):
        super().__init__()
        self.input_size = input_size
        self.output_size = output_size

        self.lin1 = nn.Linear(input_size, hidden_size[0])
        self.relu = nn.Sigmoid()
        self.lin2 = nn.Linear(hidden_size[0], hidden_size[1])
        self.relu = nn.Sigmoid()
        self.lin3 = nn.Linear(hidden_size[1], output_size)

    def forward(self, x):
        out = self.lin1(x)
        out = self.relu(out)
        out = self.lin2(out)
        out = self.relu(out)
        out = self.lin3(out)
        return out
    
        
def loss_function(Y,Y_pred):            
    Y_mask = torch.logical_not(torch.isnan(Y))      
    sq_error = torch.pow(torch.subtract(Y[Y_mask],Y_pred[Y_mask]), 2)
    mse_loss = sq_error.mean()
    return mse_loss

        
def train(X,Y):
    Y_pred = model(X)
    Y_pred = torch.squeeze(Y_pred)
    loss = loss_function(Y,Y_pred)
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()    
    return loss


@torch.no_grad()
def evaluation(x,y):
    
    Y_pred = model(x) 
    Y_pred = torch.squeeze(Y_pred)
    y = y.detach().numpy()
    Y_pred = Y_pred.detach().numpy()
    Y_mask = (~np.isnan(y))    
            
    cor = []
    mse = []
    for i in range(1):
        cor.append(np.corrcoef(y[Y_mask],Y_pred[Y_mask])[0,1])
        mse.append(np.mean((y[Y_mask]-Y_pred[Y_mask])**2))
    Cor = np.mean(cor)
    MSE = np.mean(mse)
    return float(Cor), float(MSE)


Data = CustomDataset(root= "Raw_data")

input_size = Data.X.shape[1]
hidden_size = [500,100]
output_size = 1
batch_size = 16
test_split = .2
learning_rate = 1e-4
shuffle_dataset = True

       
model = MLP(input_size, hidden_size, output_size)
  
optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)



dataset_size = Data.X.shape[0]
indices = list(range(dataset_size))
split = int(np.floor(test_split * dataset_size))

for j in (range(0,1)):
    
    if shuffle_dataset :
        np.random.shuffle(indices)
    train_indices, test_indices = indices[split:], indices[:split]
            
    train_dataset = torch.utils.data.Subset(Data, train_indices)
    test_dataset = torch.utils.data.Subset(Data, test_indices)
        
    train_loader = DataLoader(train_dataset, batch_size=batch_size)
    test_loader = DataLoader(test_dataset, batch_size=batch_size)                                               
    
    
    num_epochs = 20
    for epoch in range(num_epochs):
        print(f'Starting epoch {epoch+1}')
        
        for batch_idx, data in enumerate(train_loader):
            #print(f'batch {batch_idx}, shape {data["Input"].shape}')
            
            X = data["Input"]
            Y = data["Output"]
            X, Y = X.float(), Y.float()
            
            loss = train(X,Y)
        

        cor_train, mse_train = evaluation(train_dataset.dataset.X[train_dataset.indices],train_dataset.dataset.Y[train_dataset.indices])
        cor_test, mse_test = evaluation(test_dataset.dataset.X[test_dataset.indices],test_dataset.dataset.Y[test_dataset.indices])
        
        print(f'Epoch: {epoch:03d}, cor_train: {cor_train:.2f}, mse_train: {mse_train:.2f}, ' 
              f'cor_test: {cor_test:.2f}, mse_test: {mse_test:.2f}')
                      

        







  





          
      
#             # Print statistics
#             current_loss += loss_all_drug.item()
#             if batch_idx % 10 == 0:
#                 print('Loss after mini-batch %5d: %.3f' %
#                       (batch_idx + 1, current_loss / 500))
#                 current_loss = 0.0

# # Process is complete.
# print('Training process has finished.')
  



 


    

  
  

