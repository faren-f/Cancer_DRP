rm(list = ls())
library(e1071)
require(caTools)
library(ROCR)
library(pROC)

setwd("~/Desktop/Cancer_DRP/R/Single_Drug/")

## Read data
GE = readRDS("Raw_data/expresion_matrix.rds")
sen = readRDS("Raw_data/sensitivity_matrix.rds")

# Remove cell lines with Na from sensitivity matrix and GE matrix
i= 325
Y = sen[,i]
Y = Y[!is.na(sen[,i])]
X = GE[!is.na(sen[,i]),] # remove cell lines that are "NA" For each drug   

#External feature selection
#1) 
#Corr = cor(X,Y)
#high_corr = order(Corr,decreasing = TRUE)
#X = X[,high_corr[1:200]]

#2)
#ind_Corr = which(abs(Corr)> 0.2)
#X = X[,ind_Corr]

# GE normalization
#X = scale(X)
#X = (sel_X-min(sel_X))/(max(sel_X)-min(sel_X))


# sensitivity normalization
#Y = (Y-min(Y))/(max(Y)-min(Y))

#Binarizing drug sensitivity into sensitive or resistance
#1)
hist(Y,30)
thr = 0.6
abline(v =thr,col = "red")
Y = as.factor(ifelse(Y>thr,1,2))
table(Y)


## Classifier
## Hyperparameters
kernelParam = 0.08
Rep = 1
# From = 200
# To = 200
# Step = 1

AUC_Final = rep(0,Rep)
# for (nFeat in seq(From,To,Step)) {
AUC = rep(0,Rep)
ACC = rep(0,Rep)  
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
  
  #Undersample the majarity class
  # n1 = table(y_train)[1]
  # n2 = table(y_train)[2]
  # majorityClass = ifelse(n1>n2,1,2)
  # n_resample = abs(n1-n2)
  # 
  # ind_mC = which(y_train == levels(y_train)[majorityClass])
  # ind_mC_extra = sample(ind_mC, n_resample, replace = FALSE)
  # 
  # X_train_resampled = X_train[-ind_mC_extra,]
  # y_train_resampled = y_train[-ind_mC_extra]
  
  #X_train_resampled = X_train
  #y_train_resampled = y_train
  
  ## SVM with internal Feature selection
  #model = svm(y = y_train_resampled, x = X_train_resampled[,featureSet],kernel = "radial", gamma = kernelParam, scale = T)
  #Y_pred = predict(model, newdata = X_test[,featureSet])
  
  ## SVM without Feature selection
  model = svm(y = Ytrain, x = Xtrain, kernel = "radial", gamma = kernelParam, scale = F)
  Y_pred = predict(model, newdata = Xtest)
  Y_pred  = as.numeric(Y_pred)
  Ytest = as.numeric(Ytest)
  
  #Computing accuracy
  AC = sum((Y_pred==1 & Ytest==1) | (Y_pred==2 & Ytest==2))
  ACC[i] = AC/length(Ytest)*100

  pred = prediction(Y_pred, Ytest)
  auc = performance(pred, measure = "auc")
  AUC[i] = as.numeric(auc@y.values)
  
  #ROC curve (box plot)
  perf = performance(pred, "tpr", "fpr")
  plot(perf, avg="threshold", spread.estimate="boxplot")

  #ROC curve 
  perf = performance(pred, "tpr", "fpr")
  plot(perf, avg= "threshold", colorize=TRUE, lwd= 3,
       main= "With ROCR you can produce standard plots\nlike ROC curves ...")
  plot(perf, lty=3, col="grey78", add=TRUE)
  
  #Precision & Recall
  perf <- performance(pred, "prec", "rec")
  plot(perf, avg= "threshold", colorize=TRUE,
       lwd= 3, main= "... Precision/Recall graphs ...")
  plot(perf, lty=3, col="grey78", add=TRUE)
  
  #Sensitivity & Specificity
  perf <- performance(pred, "sens", "spec")
  plot(perf, avg= "threshold", colorize=TRUE,
       lwd= 3, main="... Sensitivity/Specificity plots ...")
  plot(perf, lty=3, col="grey78", add=TRUE)
  
  
  perf <- performance(pred, "lift", "rpp")
  plot(perf, avg= "threshold", colorize=TRUE,
       lwd= 3, main= "... and Lift charts.")
  plot(perf, lty=3, col="grey78", add=TRUE)
  
  
}

