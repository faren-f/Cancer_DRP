#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 23:01:30 2022

@author: faren
"""

import torch

def nanmean(x):
    mask = torch.logical_not(torch.isnan(x))
    len_each_col = torch.nansum(mask,dim=0) 
    mean = (torch.nansum(x,dim=0)/len_each_col)
    return(mean)