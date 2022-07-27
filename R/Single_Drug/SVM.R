rm(list = ls())
library(e1071)
require(caTools)
library(ROCR)

setwd("~/Desktop/Cancer_DRP/R/Single_Drug/")

## Read data
GE = readRDS("Data/Processed_Data/expresion_matrix.rds")
sen = readRDS("Data/Processed_Data/sensitivity_matrix.rds")

# Remove cell lines with Na from sensitivity matrix and GE matrix
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
#X_norm = scale(X)
#X_norm = (sel_X-min(sel_X))/(max(sel_X)-min(sel_X))
X_norm = X

# sensitivity normalization
Y_norm = (Y-min(Y))/(max(Y)-min(Y))

#Binarizing drug sensitivity into sensitive or resistance
#1)
# thr = 1.15
# hist(Y,30)
# abline(v =thr,col = "red")
# Y = ifelse(Y>thr_target,1,0)
# Y = as.factor(Y+1)

#2)
median_sen_mat = median(Y)
Y = as.factor(ifelse (Y>median_sen_mat,1,2))


## Classifier
## Hyperparameters
kernelParam = 0.08
Rep = 10
Result_AUC_Final = c()
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
  #model = svm(y = as.factor(y_train_resampled), x = X_train_resampled[,featureSet],kernel = "radial", gamma = kernelParam, scale = T)
  #y_hat_i = predict(model, newdata = rbind(X_test[,featureSet]))
  
  ## SVM without Feature selection
  model = svm(y = as.factor(Ytrain), x = Xtrain, kernel = "radial", gamma = kernelParam, scale = F)
  Y_pred = predict(model, newdata = rbind(X_test))
  
  #AC = sum((as.numeric(prediction)==1 & as.numeric(Ytest)==1) |
  #(as.numeric(prediction)==2 & as.numeric(Ytest)==2))
  #ACC = AC/length(Ytest)*100

  pred = prediction(y_hat, y)
  auc = performance(pred, measure = "auc")
  
  #auc = performance(pred, "tpr", "fpr")
  
  Result_AUC = c(Result_AUC, as.numeric(auc@y.values))
  #plot(Result_AUC,
  #avg= "threshold", colorize=TRUE, lwd= 3,
  #main= "With ROCR you can produce standard plots\nlike ROC curves ...")

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



