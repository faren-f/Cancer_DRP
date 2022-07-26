#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 23 7:45:38 2022

@author: Faren
"""

import numpy as np
import pandas as pd
import random
from sklearn.preprocessing import StandardScaler

import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor

# load the dataset
GE = pd.read_csv("Raw_Data/expresion_matrix.csv", header = 0, index_col=0, sep = ',' )
sen = pd.read_csv("Raw_Data/sensitivity_matrix.csv", header= None, sep = ',')

feature_list = list(GE.columns)

GE = np.array(GE)
sen = np.array(sen)
sen = sen[:,324]

Y = sen[np.where(~np.isnan(sen))].copy()
X = GE[np.where(~np.isnan(sen))[0],:].copy()
Y = np.expand_dims(Y, axis = 1)

cor_features  =   [] 
for i in range(X.shape[1]):
    cor_features.append(np.corrcoef(np.transpose(X[:,i]),np.transpose(Y))[0,1])
    
cor_features = np.abs(cor_features)
sorted_indices = np.argsort(cor_features)
reverse_sorted_indices = sorted_indices[::-1]
X = X[:,reverse_sorted_indices[0:500]]

# Normalization
ss = StandardScaler()
X = ss.fit_transform(X)
#Y = ss.fit_transform(Y)

# Instantiate model 
model = RandomForestRegressor(n_estimators = 200, random_state = 20, 
                              max_samples = 300, min_samples_split =10,
                              min_samples_leaf = 10, max_features = 200,
                              n_jobs=-1)


train_split = .8
a = range(0,X.shape[0])
num_repeat = 10
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
  
    #from sklearn.model_selection import train_test_split
    #train_features, test_features, train_labels, test_labels = train_test_split(X, Y, test_size = 0.25)

    
    # Train the model on training data
    model.fit(X_train, Y_train)

    # Use the forest's predict method on the test data
    Y_pred = model.predict(X_test)
    # Calculate the absolute errors
    Y_test = np.squeeze(Y_test)
    MSE.append(np.mean(np.square(Y_pred - Y_test)))
    Cor.append(np.corrcoef(Y_pred,Y_test)[0,1])
  
    
print(f'Mean Square Error: {np.mean(MSE)}, Cor: {np.mean(Cor)}')






# obtain n_estimators
# n_estimators = [1, 2, 4, 8, 16, 32, 64, 100, 200]
# train_results = []
# test_results = []
# for estimator in n_estimators:
#    rf = RandomForestClassifier(n_estimators=estimator, n_jobs=-1)
#    rf.fit(x_train, y_train)
#    train_pred = rf.predict(x_train)
#    false_positive_rate, true_positive_rate, thresholds = roc_curve(y_train, train_pred)
#    roc_auc = auc(false_positive_rate, true_positive_rate)
#    train_results.append(roc_auc)
#    y_pred = rf.predict(x_test)
#    false_positive_rate, true_positive_rate, thresholds = roc_curve(y_test, y_pred)
#    roc_auc = auc(false_positive_rate, true_positive_rate)
#    test_results.append(roc_auc)
# from matplotlib.legend_handler import HandlerLine2D
# line1, = plt.plot(n_estimators, train_results, ‘b’, label=”Train AUC”)
# line2, = plt.plot(n_estimators, test_results, ‘r’, label=”Test AUC”)
# plt.legend(handler_map={line1: HandlerLine2D(numpoints=2)})
# plt.ylabel(‘AUC score’)
# plt.xlabel(‘n_estimators’)
# plt.show()




# ## obtain max_depth

# max_depths = np.linspace(1, 32, 32, endpoint=True)
# train_results = []
# test_results = []
# for max_depth in max_depths:
#    rf = RandomForestClassifier(max_depth=max_depth, n_jobs=-1)
#    rf.fit(x_train, y_train)
#    train_pred = rf.predict(x_train)
#    false_positive_rate, true_positive_rate, thresholds = roc_curve(y_train, train_pred)
#    roc_auc = auc(false_positive_rate, true_positive_rate)
#    train_results.append(roc_auc)
#    y_pred = rf.predict(x_test)
#    false_positive_rate, true_positive_rate, thresholds = roc_curve(y_test, y_pred)
#    roc_auc = auc(false_positive_rate, true_positive_rate)
#    test_results.append(roc_auc)
# from matplotlib.legend_handler import HandlerLine2D
# line1, = plt.plot(max_depths, train_results, ‘b’, label=”Train AUC”)
# line2, = plt.plot(max_depths, test_results, ‘r’, label=”Test AUC”)
# plt.legend(handler_map={line1: HandlerLine2D(numpoints=2)})
# plt.ylabel(‘AUC score’)
# plt.xlabel(‘Tree depth’)
# plt.show()







