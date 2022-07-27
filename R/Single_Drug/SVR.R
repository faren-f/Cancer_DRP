rm(list = ls())
library(e1071)
require(caTools)

setwd("~/Desktop/Cancer_DRP/R/Single_Drug/")

## Read data
GE = readRDS("Data/Processed_Data/expresion_matrix.rds")
sen = readRDS("Data/Processed_Data/sensitivity_matrix.rds")

Y = sen[,i]
Y = Y[!is.na(sen[,i])]
X = GE[!is.na(sen[,i]),] # remove cell lines that are "NA" For each drug   

## External feature selection
#1) 
#Corr = cor(X,Y)
#high_corr = order(Corr,decreasing = TRUE)
#X = X[,high_corr[1:200]]

#2)
#ind_Corr = which(abs(Corr)> 0.2)
#X = X[,ind_Corr]

# GE normalization
#X_norm = scale(X)
#X_norm = (sel_X-min(sel_X))/(max(sel_X)-min(sel_X))
X_norm = X

# sensitivity normalization
Y_norm = (Y-min(Y))/(max(Y)-min(Y))


## Regression
## Hyperparameters
kernelParam = 0.08
Rep = 10
# From = 200
# To = 200
# Step = 1

Final_Cor = rep(0,Rep)
MSE = rep(0,Rep)
# for (nFeat in seq(From,To,Step)) {
#   Result_AUC = c()

for (j in 1:Rep){
  print(j)
  i=325
  
  ## Split data into train & test
  sample = sample.split(Y_norm, SplitRatio = .8)
  
  Xtrain = subset(X_norm, sample == TRUE)
  Xtest  = subset(X_norm, sample == FALSE)
  Ytrain = subset(Y_norm, sample == TRUE)
  Ytest  = subset(Y_norm, sample == FALSE)
  
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
  #Y_pred = predict(model, newdata = rbind(Xtest[,featureSet]))

  ## SVM without Feature selection
  model = svm(y = Ytrain, x = Xtrain, kernel = "radial", gamma = kernelParam, scale = F)
  Y_pred = predict(model, newdata = rbind(Xtest))
  
#}
  pred = prediction(Y_pred, Ytest)
  auc = performance(pred, measure = "auc")


  Final_Cor[j] = cor(Ytest,y_hat)
  MSE[j] = mean((Ytest-y_hat)^2)

  print(Final_Cor)[j]
  print(MSE)[j]
}
print(Final_Cor)
#print(MSE)
print(mean(Final_Cor))
#print(sd(Final_Cor))
#print(dim(X_norm))
#print(mean(MSE))
#plot(Final_Cor,MSE)



