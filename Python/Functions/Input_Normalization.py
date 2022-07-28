#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 23:01:30 2022

@author: faren

Discription: in Input_Normalization class we first normalize train data 
and then using mean and std of train data normalize test data. 
"""

import numpy as np

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
