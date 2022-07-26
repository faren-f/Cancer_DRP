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
    
@torch.no_grad()
def evaluation(x,y):
    Y_hat = model(x) 
    y = y.detach().numpy()
    Y_hat = Y_hat.detach().numpy()
    y = np.squeeze(y)
    Y_hat = np.squeeze(Y_hat)
    Cor  = np.corrcoef(y,Y_hat)[0,1]
    MSE = np.mean((y-Y_hat)**2)
    return float(Cor), float(MSE)


Data = CustomDataset(root= "Raw_data")
input_size = Data.X.shape[1]
hidden_size = [100,10]
output_size = 1
batch_size = 5
test_split = .2
learning_rate = 1e-4
shuffle_dataset = True

       
model = MLP(input_size, hidden_size, output_size)
  
loss_criterion = nn.MSELoss()
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
        trainLoss = 0
        samples = 0
        model.train()
 
        for batch_idx, data in enumerate(train_loader):
            
            X = data["Input"]
            Y = data["Output"]
            X, Y = X.float(), Y.float()
            
            Y_hat = model(X)
            loss = loss_criterion(Y, Y_hat)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step() 
            
            trainLoss += loss.item() * Y.size(0)
            samples += Y.size(0)
        trainTemplate = "epoch: {} train loss: {:.3f}"
        print(trainTemplate.format(epoch + 1, (trainLoss / samples)))
        
        
        cor_train, mse_train = evaluation(train_dataset.dataset.X[train_dataset.indices],train_dataset.dataset.Y[train_dataset.indices])
        cor_test, mse_test = evaluation(test_dataset.dataset.X[test_dataset.indices],test_dataset.dataset.Y[test_dataset.indices])
    
        print(f'Epoch: {epoch:03d}, cor_train: {cor_train:.2f}, mse_train: {mse_train:.2f}, ' 
              f'cor_test: {cor_test:.2f}, mse_test: {mse_test:.2f}')
                      

      
    
    
    
    
    
    
    
    
    
    
    
    
    
        # print(f'Epoch: {epoch:03d}, cor_train: {cor_train:.2f}, mse_train: {mse_train:.2f}, ' 
        #       f'cor_test: {cor_test:.2f}, mse_test: {mse_test:.2f}')
     
        
        
        # testLoss = 0
        # samples = 0
        # model.eval()
        # with torch.no_grad():
        #     for batch_idx, data in enumerate(test_loader):
        #         X = data["Input"]
        #         Y = data["Output"]
        #         X, Y = X.float(), Y.float()
        #         Y_hat = model(X)
        #         loss = loss_criterion(Y, Y_hat)
            
        #         testLoss += loss.item() * Y.size(0)
        #         samples += Y.size(0)
            
        #     testTemplate = "epoch: {} test loss: {:.3f}"
        #     print(testTemplate.format(epoch + 1, (testLoss / samples)))	


        
        
        
        
        
        
        
        
        
                     
