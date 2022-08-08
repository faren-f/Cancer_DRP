rm(list=ls())

library(caTools)
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
source ("RandomForest_Func.R")
source ("ENet_Func.R")
source("MLP_Func.R")
#Read data--------------------------------------------------
GE = readRDS("Processed_Data/Step1/expresion_matrix.rds")
sen = readRDS("Processed_Data/Step1/sensitivity_matrix.rds")
dim(GE)
dim(sen)

#loop across drugs--------------------------------------
Rep = 2

#mean_mse = rep(0,ncol(sen))
#sd_mse = rep(0,ncol(sen))
#mean_corr = rep(0,ncol(sen))
#sd_corr = rep(0,ncol(sen))
Results = data.frame()
for (i in 1432:1432){
  #i=325
  print(i)
  X = GE[!is.na(sen[,i]),]
  dim(X)
  y = sen[!is.na(sen[,i]),i]
  length(y)
  
  Corr = cor(X,y)
  order_corr = order(Corr,decreasing = TRUE)
  X = X[,order_corr[1:1000]]
  
  
  #loop across repeats
  mse_RF = rep(0,Rep)
  corr_RF = rep(0,Rep)
  mse_ENet = rep(0,Rep)
  corr_ENet = rep(0,Rep)
  mse_MLP = rep(0,Rep)
  corr_MLP = rep(0,Rep)
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
    
    y_pred_RF = RandomForest(ytrain = ytrain ,Xtrain = Xtrain,
                             Xtest = Xtest,ntree,mtry)
    
    y_pred_ENet = ElasticNet(ytrain = ytrain ,Xtrain = Xtrain,
                             Xtest = Xtest)
    
    y_pred_MLP = MLP(ytrain = ytrain ,Xtrain = Xtrain,
                     Xtest = Xtest)
    
    
    # y_pred re-normalization
    y_pred_RF = (y_pred_RF*STD_y)+Mean_y
    y_pred_ENet = (y_pred_ENet*STD_y)+Mean_y
    y_pred_MLP = (y_pred_MLP*STD_y)+Mean_y
    
    
    # Evaluation
    mse_RF[j] = mean((ytest-y_pred_RF)^2)
    corr_RF[j] = cor(ytest,y_pred_RF)
    
    mse_ENet[j] = mean((ytest-y_pred_ENet)^2)
    corr_ENet[j] = cor(ytest,y_pred_ENet)
    
    mse_MLP[j] = mean((ytest-y_pred_MLP)^2)
    corr_MLP[j] = cor(ytest,y_pred_MLP)
    
  }
  Results = rbind(Results, data.frame(mean_mse_RF = mean(mse_RF),
                                 sd_mse_RF = sd(mse_RF),
                                 mean_corr_RF = mean(corr_RF),
                                 sd_corr_RF = sd(corr_RF),
                                 mean_mse_ENet = mean(mse_ENet),
                                 sd_mse_ENet = sd(mse_ENet),
                                 mean_corr_ENet = mean(corr_ENet),
                                 sd_corr_ENet = sd(corr_ENet),
                                 mean_mse_MLP = mean(mse_MLP),
                                 sd_mse_MLP = sd(mse_MLP),
                                 mean_corr_MLP = mean(corr_MLP),
                                 sd_corr_MLP = sd(corr_MLP)))
  
}

Results = t(Results)
#colnames(Results) = paste0("drug_",c(1:2))
