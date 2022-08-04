rm(list=ls())

library(caTools)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
source ("RandomForest_Func.R")

#Read data--------------------------------------------------
GE = readRDS("Processed_Data/Step1/expresion_matrix.rds")
sen = readRDS("Processed_Data/Step1/sensitivity_matrix.rds")
dim(GE)
dim(sen)

#loop across drugs--------------------------------------
Rep = 5

mean_mse = rep(0,ncol(sen))
sd_mse = rep(0,ncol(sen))
mean_corr = rep(0,ncol(sen))
sd_corr = rep(0,ncol(sen))
for (i in 1:ncol(sen)){
  print(i)
  #i=325
  X = GE[!is.na(sen[,i]),]
  dim(X)
  y = sen[!is.na(sen[,i]),i]
  length(y)
  
  Corr = cor(X,y)
  high_corr = order(Corr,decreasing = TRUE)
  X = X[,high_corr[1:200]]
  
  
  #loop across repeats
  mse = rep(0,Rep)
  corr = rep(0,Rep)
  for(j in 1:Rep){
    
    sample = sample.split(y, SplitRatio = .8)
    
    Xtr_val = subset(X, sample == TRUE)
    Xtest  = subset(X, sample == FALSE)
    ytr_val = subset(y, sample == TRUE)
    ytest  = subset(y, sample == FALSE)
    
    #sample = sample.split(ytr_val, SplitRatio = .8)
    
    #Xtrain = subset(Xtr_val, sample == TRUE)
    #Xval  = subset(Xtr_val, sample == FALSE)
    #ytrain = subset(ytr_val, sample == TRUE)
    #yval  = subset(ytr_val, sample == FALSE)
    
    Xtrain = Xtr_val
    ytrain = ytr_val
    
    # Normalization
    # Xtrain normalization
    Mean_X = apply(Xtrain,2,mean)
    STD_X = apply(Xtrain,2,sd)
    Xtrain = (Xtrain-Mean_X)/STD_X

    # Xtest normalization
    Xtest = (Xtest-Mean_X)/STD_X

    # Ytrain normalization
    Mean_y = mean(ytrain)
    STD_y = sd(ytrain)
    ytrain_norm = (ytrain-Mean_y)/STD_y

    
    # Models
    ntree = 200
    mtry = 100
    y_pred = RandomForest(ytrain = ytrain ,Xtrain = Xtrain,
                      Xtest = Xtest,ntree,mtry)
    
    # y_pred re-normalization
    y_pred = (y_pred*STD_y)+Mean_y

    # Evaluation
    mse[j] = mean((ytest-y_pred)^2)
    corr[j] = cor(ytest,y_pred)

    print(mse)[j]
    print(corr)[j]
    
  }
  mean_mse[i] = mean(mse)
  sd_mse = sd(mse)
  mean_corr = mean(corr)
  sd_corr = sd(corr)
}
