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
                         repeats = 5,
                         #search = "random",
                         verboseIter = FALSE)
  
  # tune = expand.grid(alpha = seq(.05, 2, length = 20), 
  #                     lambda = round(exp(seq(-15,-5,by = 0.1)), 7))
  tune = expand.grid(alpha = seq(.05, 2, length = 20),
                     lambda = seq(.000001,0.0001,.000001))
  
  #tune = expand.grid(alpha = seq(.05, 1, length = 15),lambda = seq(.000001,.0001,.000001))

  model = caret::train(ytrain ~., data = train_data,
                         method = "glmnet",
                         metric="RMSE",
                         allowParallel = TRUE,
                         tuneGrid = tune,
                         trControl = control)
  
  y_pred = predict(model,Xtest)
  #plot(model$results$RMSE)
  #model$bestTune
  return(y_pred)
  
}


