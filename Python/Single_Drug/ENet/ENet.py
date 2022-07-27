#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 24 12:45:38 2022

@author: Faren
"""

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedKFold
from sklearn.linear_model import ElasticNet,ElasticNetCV
from sklearn.model_selection import train_test_split
from sklearn import metrics

# load the dataset
GE = pd.read_csv("Raw_data/expresion_matrix.csv", header = 0, index_col=0, sep = ',' )
sen = pd.read_csv("Raw_data/sensitivity_matrix.csv", header= None, sep = ',')

feature_list = list(GE.columns)

GE = np.array(GE)
sen = np.array(sen)
sen = sen[:,10]

Y = sen[np.where(~np.isnan(sen))].copy()
X = GE[np.where(~np.isnan(sen))[0],:].copy()
Y = np.expand_dims(Y, axis = 1)

# Normalization
ss = StandardScaler()
X = ss.fit_transform(X)
Y = ss.fit_transform(Y)

# define model
ratio = 0.06
alpha = 1.0
model = ElasticNet(alpha= alpha, l1_ratio= ratio)
#intercept = 0.035

num_repeat = 100
Cor = []
MSE = []
Score = []
for i in range(num_repeat):
    
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y,test_size = 0.1)
    
    # define model evaluation method
    #cv = RepeatedKFold(n_splits=10, n_repeats=3)
    # define model
    #ratios = np.arange(0, 1, 0.01)
    #alphas = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.0, 1.0, 10.0, 100.0]
    #model = ElasticNetCV(l1_ratio=ratios, alphas=alphas, cv=cv, n_jobs=-1)
    
    
    #fit model
    model.fit(X_train, Y_train)
    
    # summarize chosen configuration
    #print(model.alpha_)
    #print(model.intercept_)
    #print('alpha: %f' % model.alpha_)
    #print('l1_ratio_: %f' % model.l1_ratio_)
    

    # make a prediction
    Y_pred = model.predict(X_test)
    Score.append(model.score(X_test, Y_test))
    MSE.append(metrics.mean_squared_error(Y_test, Y_pred))
    Y_test = np.squeeze(Y_test)
    Cor.append(np.corrcoef(Y_test, Y_pred)[0,1])

    # summarize prediction
print(f'R2: {np.mean(Score)}, MSE: {np.mean(MSE)}, Cor: {np.mean(Cor)}')

