#                  Created on Wed Aug 4 12:04 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function receives train and test data and hyper parameters 
# of the model to compute output prediction using linear regression model 
# regularized by L2 norm popular as Ridge regression model

library(glmnet)
library(caret)


Ridge = function(ytrain,Xtrain,Xtest){
  
  train_data = cbind(Xtrain,ytrain)
  control = trainControl(method = "repeatedcv",
                         number = 5,
                         repeats = 5,
                         #search = "random",
                         verboseIter = TRUE)
  
  tune = expand.grid(alpha = 0, lambda = c(0.001,0.01,0.1,1,10))
  model = caret::train(ytrain ~., data = train_data,
                       method = "glmnet",
                       metric="RMSE",
                       allowParallel = TRUE,
                       tuneGrid = tune,
                       trControl = control)
  y_pred = predict(model,Xtest)
  
  return(y_pred)
  
}




