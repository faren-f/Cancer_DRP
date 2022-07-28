#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 17:21:42 2022

@author: faren
"""


import torch

def loss_function(Y,Y_hat):            
    Y_mask = torch.logical_not(torch.isnan(Y))      
    sq_error = torch.pow(torch.subtract(Y[Y_mask],Y_hat[Y_mask]), 2)
    mse_loss = sq_error.mean()
    return mse_loss
