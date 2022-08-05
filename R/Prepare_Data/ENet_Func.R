#                  Created on Wed Aug 4 12:04 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function receives train and test data and hyper parameters 
# of the model to compute output prediction using linear regression model 
# regularized by L1 and L2 norms popular as Elastic Net model

library(glmnet)
library(caret)

ElasticNet = function(ytrain,Xtrain,Xtest){
  
  train_data = cbind(Xtrain,ytrain)
  control = trainControl(method = "repeatedcv",
                         number = 5,
                         repeats = 1,
                         #search = "random",
                         verboseIter = TRUE)
  
  tune = expand.grid(alpha = seq(.05, 1, length = 15),
                     lambda = seq(0.001,0.1,by = 0.01))
  model = train(ytrain ~., data = train_data,
                         method = "glmnet",
                         metric="RMSE",
                         allowParallel = TRUE,
                         tuneGrid = tune,
                         trControl = control)
  
  y_pred = predict(model,Xtest)
  
  return(y_pred)
  
}


