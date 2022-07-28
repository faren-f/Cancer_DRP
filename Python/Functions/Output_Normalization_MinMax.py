#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 17:17:54 2022

@author: faren

Discription: in output_normalization function we normalize train output data 
and after training the model, using output_renormalization function 
we renormalize it. 
"""

import torch

class Output_Normalization():
    
    def __init__ (self,Ytr):
        self.Ytr = Ytr
        self.Max = torch.max(Ytr[torch.logical_not(torch.isnan(Ytr))])
        self.Min = torch.min(Ytr[torch.logical_not(torch.isnan(Ytr))])

    def normalization(self):
        
        Ytr_norm = torch.div(torch.sub(self.Ytr,self.Min), torch.sub(self.Max,self.Min))
        return Ytr_norm

    def renormalization(self):
        
        Y_tr_denorm = torch.add(torch.mul(self.Ytr_norm, (self.max-self.min)), self.min) 
        return(Y_tr_denorm)


