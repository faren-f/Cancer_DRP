#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 16:52:52 2022

@author: faren

Discription: in output_normalization function we normalize train output data 
and after training the model, using output_renormalization function 
we renormalize it. 
"""


import torch


class Output_Normalization():
    def __init__ (self,Ytr):
        
        self.Ytr = Ytr
        self.Mean = torch.mean(Ytr[torch.logical_not(torch.isnan(Ytr))])
        self.STD = torch.std(Ytr[torch.logical_not(torch.isnan(Ytr))])
        
    def normalization(self):
        
        Ytr_norm = torch.div(torch.sub(self.Ytr,self.Mean),self.STD)
        return (Ytr_norm)
    
    def renormalization(self):
        
        Ytr_denorm = torch.add(torch.mul(self.Ytr_norm, self.STD), self.Mean) 
        return (Ytr_denorm)








