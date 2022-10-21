#                  Created on Wed Oct. 9 17:20 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function receives train and test data and hyper parameters 
# of the model to compute output prediction using linear regression model 
# regularized by L2 norm popular as Ridge regression model

library(glmnet)
library(caret)

LogisticRegresion = function(ytrain, Xtrain, Xtest){
  
  train_data = cbind(Xtrain,ytrain)
  
  control = trainControl(method = 'repeatedcv',
                            number = 5,
                            repeats =  5)
  
  model = caret::train(ytrain ~., data = train_data,
                    method = 'glmnet',
                    trControl = control,
                    allowParallel = TRUE,
                    family = 'binomial')
  
  y_pred = predict(model,Xtest)
  Beta = as.matrix(coef(model$finalModel, model$bestTune$lambda))
  #order(Beta, decreasing = TRUE)[1:10]
  #plot(model$results$RMSE)
  #hist(ytrain)
  #model$bestTune
  #cor(ytest,y_pred)
  #cor(Xtest[,109],ytest)
  return(y_pred)
  
}


