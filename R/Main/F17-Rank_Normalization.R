#                  Created on Wed Sep 4 13:12 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function receives train and test data and use 
# rank normalization to normalize them

Xtrain = readRDS("Processed_data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")
Xtest = readRDS("Processed_data/S23/expresion_matrix_TCGA.rds")

Rank = function(Xtrain,Xtest){
  
  X_tr = matrix(0,nrow(Xtrain),ncol(Xtrain))
  rownames(X_tr) = rownames(Xtrain)
  colnames(X_tr) = colnames(Xtrain)
  for(j in 1:nrow(Xtrain)){
    X_tr[j,] = order(Xtrain[j,])/ncol(Xtrain)
  }
  Xtrain = X_tr

  X_te = matrix(0,nrow(Xtest),ncol(Xtest))
  rownames(X_te) = rownames(Xtest)
  colnames(X_te) = colnames(Xtest)
  for(k in 1:nrow(Xtest)){
    X_te[k,] = order(Xtest[k,])/ncol(Xtest)
  }
  Xtest = X_te
  
  X_Normalization = list(Xtrain,Xtest)
  return(X_Normalization)
}
  
  
