rm(list = ls())
library(e1071)
require(caTools)

setwd("~/Desktop/Cancer_DRP/R/Single_Drug/")

## Read data
GE = readRDS("Raw_data/expresion_matrix.rds")
sen = readRDS("Raw_data//sensitivity_matrix.rds")

i=325
Y = sen[,i]
Y = Y[!is.na(sen[,i])]
X = GE[!is.na(sen[,i]),] # remove cell lines that are "NA" For each drug   

## External feature selection
#1) 
Corr = cor(X,Y)
high_corr = order(Corr,decreasing = TRUE)
X = X[,high_corr[1:200]]

#2)
#ind_Corr = which(abs(Corr)> 0.2)
#X = X[,ind_Corr]

# GE normalization
#X = scale(X)
#X = (sel_X-min(sel_X))/(max(sel_X)-min(sel_X))

# sensitivity normalization
#Y = (Y-min(Y))/(max(Y)-min(Y))
#Y = scale(Y)

## Regression
## Hyperparameters
kernelParam = .000009
Rep = 20
# From = 200
# To = 200
# Step = 1

Final_Cor = rep(0,Rep)
MSE = rep(0,Rep)
# for (nFeat in seq(From,To,Step)) {

for (i in 1:Rep){
  print(i)

  ## Split data into train & test
  sample = sample.split(Y, SplitRatio = .8)
  
  Xtrain = subset(X, sample == TRUE)
  Xtest  = subset(X, sample == FALSE)
  Ytrain = subset(Y, sample == TRUE)
  Ytest  = subset(Y, sample == FALSE)
  
  ## Internal Feature selection
  #1)corr with sen
  # nFeat = 150
  #FeatureCorrelation = apply(X_train, 2, function(x){abs(cor(x,as.numeric(y_train)))})
  #featureSet = order(FeatureCorrelation, decreasing = TRUE)[1:nFeat]
  #hist(featureSet,30)
  
  ##2)p_val
  ## P_val distribution
  #  FeaturePval = apply(X_train, 2, function(x){
  #   a = x[y_train==levels(y_train)[1]]
  #   b = x[y_train==levels(y_train)[2]]
  #   p_val = t.test(a,b)$p.value
  #   return(p_val)
  # })
  #hist(FeaturePval,30)
  #featureSet = order(FeaturePval, decreasing = FALSE)[1:nFeat]
  
  
  ## SVM with internal Feature selection
  #model = svm(y = Ytrain, x = Xtrain[,featureSet],kernel = "radial", gamma = kernelParam, scale = T)
  #Y_pred = predict(model, newdata = Xtest[,featureSet])

  ## SVM without Feature selection
  model = svm(y = Ytrain, x = Xtrain, kernel = "radial", 
              gamma = kernelParam, scale = F)
  Y_pred = predict(model, newdata = Xtest)
  
  #Computing Correlation & MSE
  Final_Cor[i] = cor(Ytest,Y_pred)
  MSE[i] = mean((Ytest-Y_pred)^2)

  #print(Final_Cor)[i]
  #print(MSE)[i]
}
print(Final_Cor)
#print(MSE)
print(mean(Final_Cor))
#print(sd(Final_Cor))
#print(mean(MSE))
#plot(Final_Cor,MSE)


