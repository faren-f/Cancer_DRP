#                  Created on Wed Aug 4 12:04 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function receives train and test data and hyper parameters 
# of the model to compute output prediction using linear regression model 
# regularized by L2 norm popular as Ridge regression model

library(glmnet)
library(caret)

Ridge = function(ytrain, Xtrain, Xtest, weight = NULL){
  
  train_data = cbind(Xtrain,ytrain)
  control = trainControl(method = "repeatedcv",
                         number = 5,
                         repeats = 5,
                         verboseIter = FALSE)
  #tune = expand.grid(alpha = 0,lambda = seq(0.01,5,by = 0.01))
  #tune = expand.grid(alpha = 0,lambda = seq(0.000001,0.00001,by = 0.0000001))

  tune = expand.grid(alpha = 0,lambda = seq(0.01,5,by = 0.1))
  
  #tune = expand.grid(alpha = 0, lambda = round(exp(seq(-7,2.3,by = 0.1)), 4))
  #tune = expand.grid(alpha = 0, lambda = seq(.000001,0.0001,.000001))
  
  model = caret::train(ytrain ~., data = train_data,
                       method = "glmnet",
                       weights = NULL,
                       metric="RMSE",
                       allowParallel = TRUE,
                       tuneGrid = tune,
                       trControl = control)
  y_pred = predict(model,Xtest)
  #Beta = as.matrix(coef(model$finalModel, model$bestTune$lambda))
  #order(Beta, decreasing = TRUE)[1:10]
  #plot(model$results$RMSE)
  #hist(ytrain)
  #model$bestTune
  #cor(ytest,y_pred)
  #cor(Xtest[,109],ytest)
  return(y_pred)
  
}




