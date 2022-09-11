#                  Created on Wed Aug 20 19:34 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function receives train and test data and hyper parameters 
# of the model to compute output prediction using linear regression model 
# regularized by L1 norm popular as Lasso regression model

library(glmnet)
library(caret)

Lasso = function(ytrain,Xtrain,Xtest){
  
  control = trainControl(method = "repeatedcv",
                         number = 5,
                         repeats = 5,
                         verboseIter = FALSE)
  
  tune = expand.grid(alpha = 1, lambda = round(exp(seq(-7,2.3,by = 0.1)), 4))
  #tune = expand.grid(alpha = 1,lambda = seq(.000001,.0001,.000001))
  #tune = expand.grid(alpha = 1, lambda = c(0.001,0.01,0.1,1,10,100))
  model = caret::train(y= ytrain,
                       x = Xtrain,
                       method = "glmnet",
                       metric="RMSE",
                       allowParallel = TRUE,
                       tuneGrid = tune,
                       trControl = control)
  
  y_pred = predict(model,Xtest)
  model$bestTune
  return(y_pred)
  
}






