#                  Created on Wed Sep 4 13:20 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function receives cell line and clinical data and use 
# combat to remove batch effects between them and then use zscore normalization 
# to normalize them

library(sva)
Combat_Scale = function(Xtrain,Xtest){
  
  batches = rep(c('cell', 'clin'), times = c(nrow(Xtrain), nrow(Xtest)))
  batches = as.factor(batches)
  
  cell_clin = t(rbind(Xtrain, Xtest))
  
  cell_clin = ComBat(cell_clin, batches)
  
  cell_clin = t(cell_clin)
  cell_clin = scale(cell_clin)

  Xtrain = cell_clin[batches == 'cell', ]
  Xtest = cell_clin[batches == 'clin', ]

  X_Normalization = list(Xtrain,Xtest)
  return(X_Normalization)
  }
