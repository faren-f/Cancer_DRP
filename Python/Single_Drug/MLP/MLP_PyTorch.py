#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 22 12:45:38 2022

@author: Faren
"""
#import os
#os.chdir("Desktop/Codes/Cancer_DRP/Python/MTL_NN")

import numpy as np
import torch
from torch import nn
from torch.utils.data import DataLoader
from Create_CustomDataset import CustomDataset
    

class Input_Normalization():
    def __init__(self,xtr):
        self.xtr = xtr
        self.m = np.mean(self.xtr)
        self.s = np.std(self.xtr)

    def normalization_train(self):
        norm = (self.xtr-self.m)/self.s
        return norm

    def normalization_test(self,xte): 
        norm = (xte-self.m)/self.s
        return norm


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
        #self.apply(self._init_weights)
        self.relu = nn.Sigmoid()

        self.lin2 = nn.Linear(hidden_size[0], hidden_size[1])
        self.relu = nn.Sigmoid()
        #self.lin3 = nn.Linear(hidden_size[1], hidden_size[2])
        #self.relu = nn.ReLU()
        self.lin3 = nn.Linear(hidden_size[1], output_size)
        
    def _init_weights(self, m):
        if isinstance(m, nn.Linear):
            nn.init.normal_(m.weight,mean=0.0, std=1.0)
            m.bias.data.fill_(0)
            
            
    def forward(self, x):
        out = self.lin1(x)
        out = self.relu(out)
        out = self.lin2(out)
        out = self.relu(out)
        out = self.lin3(out)

        #out = self.relu(out)
        #out = self.lin4(out)
        # no activation and no softmax at the end
        return out
    

# Test the model 
# In test phase, we don't need to compute gradients (for memory efficiency)  
@torch.no_grad()
def evaluation(x,y):
    Y_hat = model(x)
   #Y_hat = torch.add(torch.mul(Y_hat, Norm_Y.STD), Norm_Y.Mean) # denormalization of Y_train_Norm 
    y = y.detach().numpy()
    Y_hat = Y_hat.detach().numpy()
    y = np.squeeze(y)
    Y_hat = np.squeeze(Y_hat)
    Cor  = np.corrcoef(y,Y_hat)[0,1]
    MSE = np.mean((y-Y_hat)**2)
    return float(Cor), float(MSE)



# Read data and seperate to Train and Test
Data = CustomDataset(root= "Raw_data")
input_size = Data.X.shape[1]
hidden_size = [100,10]
output_size = 1
batch_size = 5
test_split = .2
learning_rate = 1e-4
shuffle_dataset = True

       
# Initialize the MLP
model = MLP(input_size, hidden_size, output_size)
  
# loss function and Optimizer
loss_function = nn.MSELoss()
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
    
    # Norm_X = Input_Normalization(train_dataset.dataset.X[train_dataset.indices])
    # X_tr_norm = Norm_X.normalization_train()
    # train_dataset.dataset.X[train_dataset.indices] = X_tr_norm
    
    # X_te_norm = Norm_X.normalization_test(test_dataset.dataset.X[test_dataset.indices])
    # test_dataset.dataset.X[test_dataset.indices] = X_te_norm
    
    # Norm_Y = Output_Normalization(train_dataset.dataset.Y[train_dataset.indices])
    # Y_tr_norm = Norm_Y.normalization_train()
    # train_dataset.dataset.Y[train_dataset.indices] = Y_tr_norm
    
    
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
    num_epochs = 30
    for epoch in range(num_epochs):
        # Print epoch
        print(f'Starting epoch {epoch+1}')
        
        # Set current loss value
        trainLoss = 0
        samples = 0
        model.train()
        
        # Iterate over the DataLoader for training data    
        for batch_idx, data in enumerate(train_loader):
            #print(f'batch {batch_idx}, shape {data["Input"].shape}')
            
            # Get and prepare inputs
            X = data["Input"]
            Y = data["Output"]
            X, Y = X.float(), Y.float()
            
            # Perform forward pass
            Y_hat = model(X)
                
            # Compute loss
            loss = loss_function(Y, Y_hat)

            # Zero the gradients
            optimizer.zero_grad()
                
            # Perform backward pass
            loss.backward()
                
            # Perform optimization
            optimizer.step()    

            trainLoss += loss.item() * Y.size(0)
            samples += Y.size(0)
        trainTemplate = "epoch: {} train loss: {:.3f}"
        print(trainTemplate.format(epoch + 1, (trainLoss / samples)))
        
        
        #Y_tr_denorm = torch.add(torch.mul(Y_tr_norm, Norm_Y.STD), Norm_Y.Mean) # denormalization of Y_train_Norm 

        #cor_train, mse_train = evaluation(X_tr_norm,Y_tr_denorm)
        #cor_test, mse_test = evaluation(X_te_norm,test_dataset.dataset.Y[test_dataset.indices])
        
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
  



 


    

  
  