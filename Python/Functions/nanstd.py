#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 23:01:30 2022

@author: faren
"""
import torch

def nanstd(x):
    mask = torch.logical_not(torch.isnan(x))
    len_each_col = torch.nansum(mask,dim=0)
    Mean = (torch.nansum(x,dim=0)/len_each_col)
    Mean_broadcasting = Mean.repeat(x.shape[0],1)
    Subtraction = x - Mean_broadcasting
    Pow = torch.pow(Subtraction,2)
    std = (torch.nansum(Pow,dim=0)/len_each_col)
    return(std)
