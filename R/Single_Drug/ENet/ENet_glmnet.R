rm(list = ls())

require(caTools)
library(corrplot)
library(glmnet)
library(caret)

setwd("~/Desktop/Cancer_DRP/R/Single_Drug/ENet/")

## Read data
GE = readRDS("Raw_data/expresion_matrix.rds")
sen = readRDS("Raw_data/sensitivity_matrix.rds")

i=325
X = GE[!is.na(sen[,i]),]           # remove cell lines that are "NA" For each drug   
y = sen[!is.na(sen[,i]),i]

Corr = cor(X,y)
order_Corr = order(Corr, decreasing = TRUE)
X = X[,order_Corr[1:5000]]

Rep = 2
corr = rep(0,Rep)
mse = rep(0,Rep)
for (j in 1:Rep){
  print(j)
  
  ## Split data into train & test
  sample = sample.split(y, SplitRatio = .8)
  
  Xtrain = subset(X, sample == TRUE)
  Xtest  = subset(X, sample == FALSE)
  ytrain = subset(y, sample == TRUE)
  ytest  = subset(y, sample == FALSE)
  
  # Normalization
  # Xtrain normalization
  Mean_X = apply(Xtrain,2,mean)
  STD_X = apply(Xtrain,2,sd)
  Xtrain = (Xtrain-Mean_X)/STD_X

  # Xtest normalization
  Xtest = (Xtest-Mean_X)/STD_X

  #Ytrain normalization
  #Ytrain = (Ytrain-min(Ytrain))/(max(Ytrain)-min(Ytrain))
  Mean_y = mean(ytrain)
  STD_y = sd(ytrain)
  ytrain_norm = (ytrain-Mean_y)/STD_y

  train_data = cbind(Xtrain,ytrain_norm)
  
  ## Model training
  # Model Building : Elastic Net Regression
  control = trainControl(method = "repeatedcv",
                          number = 10,
                          repeats = 10,
                          #search = "random",
                          verboseIter = TRUE)
  
  tune = expand.grid(alpha = seq(.05, 1, length = 15),
                     lambda = seq(0.001,0.1,by = 0.01))
  # Training ELastic Net Regression model
  model = train(ytrain_norm ~., data = train_data,
                         method = "glmnet",
                         metric="RMSE",
                         allowParallel = TRUE,
                         tuneGrid = tune,
                         trControl = control)
  
  model$bestTune
  #coef(model$finalModel, model$bestTune$lambda)
  #alpha = 0.9           # or alpha = 1
  #lambda = 0.09
  # Model Prediction
  y_pred = predict(model, Xtest)
  
  # y_pred re-normalization
  y_pred = (y_pred*STD_y)+Mean_y
  
  mse[j] = mean((ytest - y_pred)^2)
  corr[j] = cor(ytest,y_pred, method = "pearson")
  
  print(mse)
  print(corr)
}

#print(mean(mse))
print(mean(corr))
#print(sd(corr))
