rm(list=ls())
library(mRMRe)
library(randomForest)
require(caTools)

setwd("~/Desktop/Cancer_DRP/R/Single_Drug/")

## Read data
GE = readRDS("Raw_data/expresion_matrix.rds")
sen = readRDS("Raw_data/sensitivity_matrix.rds")

i=325
Y = sen[,i]
Y = Y[!is.na(sen[,i])]
X = GE[!is.na(sen[,i]),] # remove cell lines that are "NA" For each drug   
#Corr = cor(X,Y)
#high_corr = order(Corr,decreasing = TRUE)
#X = X[,high_corr[1:200]]



###mRMR Feature selection






# Normalization
X = scale(X)
#X = (X-min(X))/(max(X)-min(X))

#Y = (Y-min(Y))/(max(Y)-min(Y))
Y = scale(Y)

## Classifier
Final_Cor = rep(0,20)
MSE = rep(0,20)

for (j in 1:20){
  print(j)

  ## Split data into train & test
  sample = sample.split(Y, SplitRatio = .8)
  
  Xtrain = subset(X, sample == TRUE)
  Xtest  = subset(X, sample == FALSE)
  Ytrain = subset(Y, sample == TRUE)
  Ytest  = subset(Y, sample == FALSE)
  
  RF = randomForest(y = Ytrain,x = Xtrain, ntree = 200,mtry = 100)
  y_hat = predict(RF, newdata=Xtest)
  
  
  Final_Cor[j] = cor(Ytest,y_hat)
  MSE[j] = mean((Ytest-y_hat)^2)
  
  print(Final_Cor)[j]
  print(MSE)[j]
}
#}
print(Final_Cor)
#print(MSE)
print(mean(Final_Cor))
#print(sd(Final_Cor))
#print(mean(MSE))
#plot(Final_Cor,MSE)



