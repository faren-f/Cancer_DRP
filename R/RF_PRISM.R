rm(list = ls())
library(randomForest)
library(caTools)

setwd("~/Desktop/Codes/Cancer_DRP/R")
## Read data
GE1 = readRDS("Data/Processed_Data/expresion_matrix.rds")
GE =  read.csv("Data/Processed_Data/expresion_matrix.csv")
sen = readRDS("Data/Processed_Data/sensitivity_matrix.rds")

#i=1429
i = 325
Y = sen[,i]
Y = Y[!is.na(sen[,i])]
X = GE[!is.na(sen[,i]),] # remove cell lines that are "NA" For each drug   

## Corr input & output
Corr = cor(X,Y)
#hist(abs(Corr))  



sample = sample.split(Y, SplitRatio = .8)

## Classifier
Final_Cor = rep(0,20)
MSE = rep(0,20)

for (j in 1:20){
  print(j)
  
  ## Split data into train & test
  Xtrain = subset(X, sample == TRUE)
  Xtest  = subset(X, sample == FALSE)
  Ytrain = subset(Y, sample == TRUE)
  Ytest  = subset(Y, sample == FALSE)
  
  RF = randomForest(y = Ytrain,x = Xtrain, ntree = 200,mtry = 100)
  Y_hat = predict(RF, newdata=Xtest)
  
  Final_Cor[j] = cor(Ytest,Y_hat)
  MSE[j] = mean((Ytest-Y_hat)^2)
  
  print(Final_Cor)[j]
  print(MSE)[j]
}
print(Final_Cor)
print(MSE)
print(mean(Final_Cor))
print(sd(Final_Cor))
print(mean(MSE))

