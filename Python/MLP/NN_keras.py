#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 23 12:45:38 2022

@author: Faren
"""

import numpy as np
import pandas as pd
import random
from sklearn.preprocessing import StandardScaler
import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
import matplotlib.pyplot as plt

# load the dataset


GE = pd.read_csv("Raw_Data/expresion_matrix.csv", header = 0, index_col=0, sep = ',' )
sen = pd.read_csv("Raw_Data/sensitivity_matrix.csv", header= None, sep = ',')


GE = np.array(GE)
sen = np.array(sen)
sen = sen[:,324]


Y = sen[np.where(~np.isnan(sen))].copy()
X = GE[np.where(~np.isnan(sen))[0],:].copy()
Y = np.expand_dims(Y, axis = 1)

# cor_features  =   [] 
# for i in range(X.shape[1]):
#     cor_features.append(np.corrcoef(np.transpose(X[:,i]),np.transpose(Y))[0,1])
    
# cor_features = np.abs(cor_features)
# sorted_indices = np.argsort(cor_features)
# reverse_sorted_indices = sorted_indices[::-1]
# X = X[:,reverse_sorted_indices[0:500]]

# Normalization
ss = StandardScaler()
X = ss.fit_transform(X)
#Y = ss.fit_transform(Y)

train_split = .8
a = range(0,X.shape[0])
num_repeat = 1
Cor = []
MSE = []
for i in range(num_repeat):
    indeces = random.sample(a,len(a))
    Tr_indeces = indeces[0:round(len(a)*train_split)]
    Te_indeces = indeces[round(len(a)*train_split):]
        
    X_train = X[Tr_indeces]
    X_test = X[Te_indeces]
    Y_train = Y[Tr_indeces]
    Y_test = Y[Te_indeces]
    
    
    input_size = X.shape[1]
    hidden_size = [100,10]
    output_size = 1
    batch_size = 10
    test_split = .2
    learning_rate = 1e-4
    
    # define the keras model
    model = Sequential()
    model.add(Dense(hidden_size[0], input_dim=input_size, 
                    activation='sigmoid'))
                    #,kernel_initializer='normal'))
    model.add(Dense(hidden_size[1], activation='sigmoid'))
    model.add(Dense(output_size, activation='linear'))
    
    
    # compile the keras model
    opt = keras.optimizers.Adam(learning_rate=0.0001)
    model.compile(loss='mse', 
                  optimizer= opt)
    
    
    # fit the keras model on the dataset
    model.fit(X_train, Y_train, epochs=20, batch_size=10)
    
    
    # evaluate the keras model
    
    Y_pred = model.predict(X_test)
    
    Cor.append(np.corrcoef(np.transpose(Y_test),np.transpose(Y_pred))[0,1])
    MSE.append(np.mean((Y_test-Y_pred)**2))
    
      
print(f'Cor: {Cor}')
print(f'MSE: {MSE}')



