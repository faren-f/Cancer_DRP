import os
#os.chdir("Desktop/Codes/Cancer_DRP/Python/MTL_NN")

import numpy as np

import torch
from torch import nn

from torch.utils.data import DataLoader
from Create_CustomDataset import CustomDataset
    

def nanmean(x):
    mask = torch.logical_not(torch.isnan(x))
    len_each_col = torch.nansum(mask,dim=0) # na is FALSE so they do not considered to find length of each column
    mean = (torch.nansum(x,dim=0)/len_each_col)
    return(mean)

def nanstd(x):
    mask = torch.logical_not(torch.isnan(x))
    len_each_col = torch.nansum(mask,dim=0)
    Mean = (torch.nansum(x,dim=0)/len_each_col)
    Mean_broadcasting = Mean.repeat(x.shape[0],1)
    Subtraction = x - Mean_broadcasting
    Pow = torch.pow(Subtraction,2)
    std = (torch.nansum(Pow,dim=0)/len_each_col)
    return(std)

 
class Input_Normalization():
    def __init__(self,xtr):
        self.xtr = xtr
        self.m = nanmean(self.xtr)
        self.s = nanstd(self.xtr)

    def normalization_train(self):
        norm = (self.xtr-self.m)/self.s
        return norm

    def normalization_test(self,xte): 
        norm = (xte-self.m)/self.s
        return norm

# Min-Max normalization(It did not have good performance)    
# def output_normalization(Ytr):
#     Max = torch.max(Ytr[torch.logical_not(torch.isnan(Ytr))])
#     Min = torch.min(Ytr[torch.logical_not(torch.isnan(Ytr))])
#     Ytr_norm = torch.div(torch.sub(Ytr,Min), torch.sub(Max,Min))
#     return Ytr_norm



class Output_Normalization():
    def __init__ (self,Ytr):
        self.Ytr = Ytr
        self.Mean = torch.mean(Ytr[torch.logical_not(torch.isnan(Ytr))])
        self.STD = torch.std(Ytr[torch.logical_not(torch.isnan(Ytr))])
    def normalization_train(self):
        norm = torch.div(torch.sub(self.Ytr,self.Mean),self.STD)
        return norm


class MLP(nn.Module):
    def __init__(self, input_size, hidden_size, output_size):
        super().__init__()
        self.input_size = input_size
        self.output_size = output_size

        self.lin1 = nn.Linear(input_size, hidden_size[0])
        self.relu = nn.ReLU()
        self.lin2 = nn.Linear(hidden_size[0], hidden_size[1])
        self.relu = nn.ReLU()
        self.lin3 = nn.Linear(hidden_size[1], output_size)

    def forward(self, x):
        out = self.lin1(x)
        out = self.relu(out)
        out = self.lin2(out)
        out = self.relu(out)
        out = self.lin3(out)
        # no activation and no softmax at the end
        return out
    
        
def loss_function(Y,Y_hat):            
    Y_mask = torch.logical_not(torch.isnan(Y))      
    sq_error = torch.pow(torch.subtract(Y[Y_mask],Y_hat[Y_mask]), 2)
    mse_loss = sq_error.mean()
    return mse_loss

        
def train(X,Y):
    # Perform forward pass
    Y_hat = model(X)
    
    # Compute loss
    loss = loss_function(Y,Y_hat)
    
    # Zero the gradients
    optimizer.zero_grad()
    
    # Perform backward pass
    loss_all_drug.backward()
    
    # Perform optimization
    optimizer.step()    
    return loss

    

# Test the model
# In test phase, we don't need to compute gradients (for memory efficiency)
@torch.no_grad()
def evaluation(x,y):
    
    Y_hat = model(x) 
    #Y_Na_zero  = y.clone().detach()          
    #Y_Na_zero = torch.where(torch.isnan(Y_Na_zero),torch.tensor(0.),y)
    Y_hat = torch.add(torch.mul(Y_hat, Norm_Y.STD), Norm_Y.Mean) # denormalization of Y_train_Norm 
    y = y.detach().numpy()
    Y_hat = Y_hat.detach().numpy()
    Y_mask = (~np.isnan(y))
    #Y_Na_zero  = Y_Na_zero.detach().numpy()
    
            
    cor = []
    mse = []
    for i in range(x.shape[0]):
        cor.append(np.corrcoef(y[Y_mask[:,i],i],Y_hat[Y_mask[:,i],i])[0,1])
        #cor.append(np.corrcoef(Y_Na_zero[i,:],Y_hat[i,:])[0,1])
        mse.append(np.mean((y[Y_mask[:,i],i]-Y_hat[Y_mask[:,i],i])**2))
    Cor = np.mean(cor)
    MSE = np.mean(mse)
    return float(Cor), float(MSE)


# Read data and seperate to Train and Test
Data = CustomDataset(root= "Raw_data")
print(Data)
input_size = Data.X.shape[1]
hidden_size = [500,100]
output_size = Data.Y.shape[1]
batch_size = 16
test_split = .2
learning_rate = 1e-4
shuffle_dataset = True

       
# Initialize the MLP
model = MLP(input_size, hidden_size, output_size)
  
# Optimizer
optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)



# Creating data indices for training and validation splits:
dataset_size = Data.X.shape[0]
indices = list(range(dataset_size))
split = int(np.floor(test_split * dataset_size))

for j in (range(0,1)):
    
    if shuffle_dataset :
        np.random.shuffle(indices)
    train_indices, test_indices = indices[split:], indices[:split]
            
    train_dataset = torch.utils.data.Subset(Data, train_indices)
    test_dataset = torch.utils.data.Subset(Data, test_indices)
    
    Norm_X = Input_Normalization(train_dataset.dataset.X[train_dataset.indices])
    X_tr_norm = Norm_X.normalization_train()
    train_dataset.dataset.X[train_dataset.indices] = X_tr_norm
    
    X_te_norm = Norm_X.normalization_test(test_dataset.dataset.X[test_dataset.indices])
    test_dataset.dataset.X[test_dataset.indices] = X_te_norm
    
    Norm_Y = Output_Normalization(train_dataset.dataset.Y[train_dataset.indices])
    Y_tr_norm = Norm_Y.normalization_train()
    train_dataset.dataset.Y[train_dataset.indices] = Y_tr_norm
    
    
    # Norm_Y = Normalization(train_dataset.dataset.Y[train_dataset.indices])
    # Y_tr_norm = Norm_Y.normalization_train()
    # train_dataset.dataset.Y[train_dataset.indices] = Y_tr_norm
    
    # Y_te_norm = Norm_Y.normalization_test(test_dataset.dataset.Y[test_dataset.indices])
    # test_dataset.dataset.Y[test_dataset.indices] = Y_te_norm
    
    
    #Method 2 for normalization from Pytorch (I am not sure about "Na" in this method)
    # from sklearn.preprocessing import StandardScaler

    # scaler = StandardScaler()
    # X_tr_norm = scaler.fit_transform(train_dataset.dataset.X[train_dataset.indices])
    # X_tr_norm = torch.tensor(X_tr_norm,dtype=torch.float32)
    # train_dataset.dataset.X[train_dataset.indices] = X_tr_norm

    # X_te_norm = scaler.transform(test_dataset.dataset.X[test_dataset.indices])
    # X_te_norm = torch.tensor(X_te_norm,dtype=torch.float32)
    # test_dataset.dataset.X[test_dataset.indices] = X_te_norm

    # Y_tr_norm = scaler.fit_transform( train_dataset.dataset.Y[train_dataset.indices] )
    # Y_tr_norm = torch.tensor(Y_tr_norm,dtype=torch.float32)
    # train_dataset.dataset.Y[train_dataset.indices] = Y_tr_norm

    # Y_te_norm = scaler.transform( test_dataset.dataset.Y[test_dataset.indices] )
    # Y_te_norm = torch.tensor(Y_te_norm,dtype=torch.float32)
    # test_dataset.dataset.Y[test_dataset.indices] = Y_te_norm

        
    train_loader = DataLoader(train_dataset, batch_size=batch_size)
    test_loader = DataLoader(test_dataset, batch_size=batch_size)                                               
    
    
    # Training loop:       
    num_epochs = 50
    for epoch in range(num_epochs):
        # Print epoch
        print(f'Starting epoch {epoch+1}')
        
        # Set current loss value
        #current_loss = 0.0
        
        # Iterate over the DataLoader for training data    
        for batch_idx, data in enumerate(train_loader):
            #print(f'batch {batch_idx}, shape {data["Input"].shape}')
            
            # Get and prepare inputs
            X = data["Input"]
            Y = data["Output"]
            X, Y = X.float(), Y.float()
            
            loss = train(X,Y)
        
        Y_tr_denorm = torch.add(torch.mul(Y_tr_norm, Norm_Y.STD), Norm_Y.Mean) # denormalization of Y_train_Norm 

        cor_train, mse_train = evaluation(X_tr_norm,Y_tr_denorm)
        cor_test, mse_test = evaluation(X_te_norm,test_dataset.dataset.Y[test_dataset.indices])
        
        # cor_train, mse_train = evaluation(X_tr_norm,train_dataset.dataset.Y[train_dataset.indices])
        # cor_test, mse_test = evaluation(X_te_norm,test_dataset.dataset.Y[test_dataset.indices])
        
        
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
  



 


    

  
  