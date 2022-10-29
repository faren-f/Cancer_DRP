#                  Created on Wed Oct. 29 17:34 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function receives train and test data to compute output 
#prediction using linear regression model 


library(glmnet)
library(caret)

LinearcRegresion = function(ytrain, Xtrain, Xtest){
  
  train_data = cbind(Xtrain,ytrain)
  
  control = trainControl(method = 'repeatedcv',
                         number = 5,
                         repeats =  5)
  
  model = caret::train(ytrain ~., data = train_data,
                       method = 'lm',
                       trControl = control,
                       allowParallel = TRUE)
  
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


