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
from sklearn.model_selection import train_test_split



## Main section

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

# Instantiate model 
model = RandomForestRegressor(n_estimators = 200, max_samples = 300, 
                              min_samples_split =10, min_samples_leaf = 10, 
                              max_features = 200, n_jobs=-1)

test_size = 0.2
num_repeat = 10
Cor = []
MSE = []
for i in range(num_repeat):
  
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size = 0.25)

    
    # Train the model on training data
    model.fit(X_train, Y_train)

    # Use the forest's predict method on the test data
    Y_pred = model.predict(X_test)
    # Calculate the absolute errors
    Y_test = np.squeeze(Y_test)
    MSE.append(np.mean(np.square(Y_pred - Y_test)))
    Cor.append(np.corrcoef(Y_pred,Y_test)[0,1])
   
print(f'Mean Square Error: {np.mean(MSE)}, Cor: {np.mean(Cor)}')




### obtain n_estimators
n_estimators = [180,190,200,210,220,230]

Cor_train = []
Cor_test = []
for estimator in n_estimators:
    model = RandomForestRegressor(n_estimators=estimator, n_jobs=-1)
    model.fit(X_train, Y_train)
    Y_train_pred = model.predict(X_train)
    #MSE.append(np.mean(np.square(Y_train_pred - Y_train)))
    Y_train = np.squeeze(Y_train)
    Cor_train.append(np.corrcoef(Y_train_pred,Y_train)[0,1])
    Y_pred = model.predict(X_test)
    Y_test = np.squeeze(Y_test)
    Cor_test.append(np.corrcoef(Y_pred,Y_test)[0,1])
  
    
#from matplotlib.legend_handler import HandlerLine2D
line1 = plt.plot(n_estimators, Cor_train, 'b', label='Cor_train')
line2 = plt.plot(n_estimators, Cor_test, 'r', label='Cor_test')
#plt.legend(handler_map={line1: HandlerLine2D(numpoints=2)})
plt.ylabel('Cor')
plt.xlabel('n_estimators')
plt.show()




### obtain max_depth

max_depths = np.linspace(1, 32, 32, endpoint=True)
Cor_train = []
Cor_test = []
for max_depth in max_depths:
    model = RandomForestRegressor(max_depth=max_depth, n_jobs=-1)
    model.fit(X_train, Y_train)
    Y_train_pred = model.predict(X_train)
    Y_train = np.squeeze(Y_train)
    Cor_train.append(np.corrcoef(Y_train_pred,Y_train)[0,1])
    Y_pred = model.predict(X_test)
    Y_test = np.squeeze(Y_test)
    Cor_test.append(np.corrcoef(Y_pred,Y_test)[0,1])

#from matplotlib.legend_handler import HandlerLine2D    
line1 = plt.plot(max_depths, Cor_train, 'b', label='Cor_train')
line2 = plt.plot(max_depths, Cor_test, 'r', label='Cor_test')
#plt.legend(handler_map={line1: HandlerLine2D(numpoints=2)})
plt.ylabel('Cor')
plt.xlabel('Tree depth')
plt.show()
    



# Import tools needed for visualization
from sklearn.tree import export_graphviz
import pydot
# Pull out one tree from the forest
tree = model.estimators_[5]

# Export the image to a dot file
export_graphviz(tree, out_file = 'tree.dot', feature_names = feature_list, rounded = True, precision = 1)
# Use dot file to create a graph
(graph, ) = pydot.graph_from_dot_file('tree.dot')
# Write graph to a png file
graph.write_png('tree.png')




# Limit depth of tree to 3 levels
model_small = RandomForestRegressor(n_estimators=10, max_depth = 3)
model_small.fit(X_train, Y_train)
# Extract the small tree
tree_small = model_small.estimators_[5]
# Save the tree as a png image
export_graphviz(tree_small, out_file = 'small_tree.dot', 
                feature_names = feature_list, rounded = True, precision = 1)
(graph, ) = pydot.graph_from_dot_file('small_tree.dot')
graph.write_png('small_tree.png');







# Get numerical feature importances
importances = list(model.feature_importances_)
# List of tuples with variable and importance
feature_importances = [(feature, round(importance, 2)) for feature, importance in zip(feature_list, importances)]
# Sort the feature importances by most important first
feature_importances = sorted(feature_importances, key = lambda x: x[1], reverse = True)
# Print out the feature and importances 
[print('Variable: {:20} Importance: {}'.format(*pair)) for pair in feature_importances];




# New random forest with only the two most important variables
model_most_important = RandomForestRegressor(n_estimators= 1000, random_state=42)
# Extract the two most important features
important_indices = [feature_list.index('temp_1'), feature_list.index('average')]
train_important = X_train[:, important_indices]
test_important = X_test[:, important_indices]
# Train the random forest
model_most_important.fit(train_important, Y_train)
# Make predictions and determine the error
predictions = model_most_important.predict(test_important)
errors = abs(predictions - Y_test)
# Display the performance metrics
print('Mean Absolute Error:', round(np.mean(errors), 2), 'degrees.')
mape = np.mean(100 * (errors / Y_test))
accuracy = 100 - mape
print('Accuracy:', round(accuracy, 2), '%.')









