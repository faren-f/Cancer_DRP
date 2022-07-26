#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 23:01:30 2022

@author: faren
"""
import torch
import numpy as np

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


# Z-score normalization
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

#Min-Max normalization   
def output_normalization(Ytr):
    Max = torch.max(Ytr[torch.logical_not(torch.isnan(Ytr))])
    Min = torch.min(Ytr[torch.logical_not(torch.isnan(Ytr))])
    Ytr_norm = torch.div(torch.sub(Ytr,Min), torch.sub(Max,Min))
    return Ytr_norm




## computing loss when we have nan in our data
def loss_function(Y,Y_hat):            
    Y_mask = torch.logical_not(torch.isnan(Y))      
    sq_error = torch.pow(torch.subtract(Y[Y_mask],Y_hat[Y_mask]), 2)
    mse_loss = sq_error.mean()
    return mse_loss


