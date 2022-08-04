
#                  Created on Wed Aug 3 12:04 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function receives train and test and Random Forest 
# hyper parameters and compute output prediction

library(randomForest)

RandomForest = function(ytrain,Xtrain,Xtest,ntree,mtry){

  RF = randomForest(y = ytrain,x = Xtrain, ntree=ntree,mtry=mtry)
  y_pred = predict(RF, newdata=Xtest)
  
  return(y_pred)
}