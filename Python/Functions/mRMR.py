
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 16:16:48 2022

@author: faren

Description: This is the mRMR function that uses sklearn package for estimating the 
distribusion of input and output using KNN for computing mutual information, and then 
we find features that have maximum relevance with output while minimum redundancy with
priviously selected features.  
"""


import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.feature_selection as skfs


# load the dataset
GE = pd.read_csv("Raw_Data/expresion_matrix.csv", header = 0, index_col=0, sep = ',' )
sen = pd.read_csv("Raw_Data/sensitivity_matrix.csv", header= None, sep = ',')


GE = np.array(GE)
sen = np.array(sen)
sen = sen[:,324]

Y = sen[np.where(~np.isnan(sen))].copy()
X = GE[np.where(~np.isnan(sen))[0],:].copy()


def mRMR(X, Y, N_feat=np.nan, alpha=1, plot= True):
    
    
    if pd.isnull(N_feat):
        N_feat = X.shape[1]
        
    
    MI_Xy = skfs.mutual_info_regression(X, Y, discrete_features='auto', 
                                                    n_neighbors=3)
    MI_XX = np.zeros([X.shape[1],X.shape[1]])
    for i in range(0,X.shape[1]):
            
            Z = X[:,i]
            #Z = np.expand_dims(Z, axis = 1)
            MI_XX[i,:] = skfs.mutual_info_regression(X, Z,
                                                discrete_features='auto',
                                                n_neighbors=3)
                 
    MI_Xy_orders = np.argsort(MI_Xy)
    MI_Xy_orders = MI_Xy_orders[::-1]
    selected_feat = MI_Xy_orders[0]
    residual_feat = MI_Xy_orders[1::]
    #residual_feat = residual_feat.tolist()
    scores = max(MI_Xy)
    alpha = 1
    while(len(residual_feat)!=0):
        
        S = []
        for j in residual_feat:
            print(j)         
            score_j = MI_Xy[j] - alpha*(np.sum(MI_XX[j, selected_feat])/selected_feat.size)
            S.append(score_j)
        S = np.array(S)
        selected_feat = np.append(selected_feat, residual_feat[np.argmax(S)])
        residual_feat = np.delete(residual_feat,np.argmax(S))
    
        scores = np.append(scores, np.max(S))
        selected_feat = selected_feat[0:N_feat]
        if plot:
            plt.plot(scores)
            
    return(selected_feat, scores)
    
        

selected_feat, scores = mRMR(X = X[0:10,0:5], Y = Y[0:10], 
                             N_feat=np.nan, alpha=1, plot= True)








