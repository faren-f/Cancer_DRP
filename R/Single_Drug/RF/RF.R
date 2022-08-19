rm(list = ls())
library(randomForest)
require(caTools)

setwd("~/Desktop/Cancer_DRP/R/Single_Drug/RF/")

## Read data
sen = readRDS("Raw_data/sensitivity_matrix.rds")
GE = readRDS("Raw_data/expresion_matrix.rds")
drug_targets = readRDS("Raw_data/drug_targets.rds")
#i = 540
i = 1020

#for(i in 1:ncol(sen)){
y = sen[,i]
y = y[!is.na(sen[,i])]
X = GE[!is.na(sen[,i]),] 

#Corr = cor(X,y)
#high_corr = order(Corr,decreasing = TRUE)
#X = X[,high_corr[1:600]]

X = scale(X)
#X = (X-min(X))/(max(X)-min(X))

y = scale(y)
#y = (y-min(y))/(max(y)-min(y))

Rep = 1
corr = rep(0,Rep)
mse = rep(0,Rep)

for (j in 1:Rep){
  print(paste0("Rep is: ",j))
  
  ## Split data into train & test
  sample = sample.split(y, SplitRatio = .8)
  
  Xtrain = subset(X, sample == TRUE)
  Xtest  = subset(X, sample == FALSE)
  ytrain = subset(y, sample == TRUE)
  ytest  = subset(y, sample == FALSE)
  
  ytrain = as.vector(ytrain)
  ytest = as.vector(ytest)
  
  model = randomForest(y = ytrain,x = Xtrain, ntree = 200,mtry = 100)
  y_pred = predict(model, newdata=Xtest)
  
  corr[j] = cor(ytest,y_pred)
  mse[j] = mean((ytest-y_pred)^2)
  
  print(corr)
  #print(mse)
}
#}
#print(corr)
#print(mse)
print(mean(corr))
#print(sd(corr))
print(mean(mse))
#plot(corr,mse)



