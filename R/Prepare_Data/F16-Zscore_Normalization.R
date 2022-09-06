#                  Created on Wed Sep 4 13:01 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function receives train and test data and use 
# zscore normalization to normalize them

ZScore = function(Xtrain,Xtest){
  
  Xtr = t(Xtrain)
  Xtr = scale(Xtr)
  Xtrain = t(Xtr)
  
  Xte = t(Xtest)
  Xte = scale(Xte)
  Xtest = t(Xte)
  
  X_Normalization = list(Xtrain, Xtest)
  
  return(X_Normalization)
}





