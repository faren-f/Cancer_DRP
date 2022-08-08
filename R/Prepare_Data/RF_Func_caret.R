
#                  Created on Wed Aug 8 12:04 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function receives train and test and hyper parameters of the model 
# to compute output prediction using Random Forest 

library(randomForest)
library(caret)
RandomForest = function(ytrain,Xtrain,Xtest){
  
  train_data = cbind(Xtrain,ytrain)
  
  control = trainControl(method = "cv",
                         number = 10,
                         search = "grid",
                         allowParallel = TRUE)
  
  
  tuneGrid = expand.grid(.mtry = seq(50, 200, by=5))
  
  maxtrees = list()
  for (ntree in c(50,70,100,150,200,250,300)) {
    model = caret::train(ytrain~., 
                          data = train_data,
                          method = "rf",
                          tuneGrid = tuneGrid,
                          trControl = control,
                          ntree = ntree)
    key = toString(ntree)
    maxtrees[[key]] = model
  }
  results_tree = resamples(maxtrees)
  
  y_pred = predict(model, Xtest)
  
  return(y_pred)
}

