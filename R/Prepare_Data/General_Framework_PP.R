rm(list=ls())

library(caTools)
library(parallel)

no_cores = detectCores()
cl = makeCluster(no_cores-2)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
source ("RandomForest_Func.R")

#Read data--------------------------------------------------
GE = readRDS("Processed_Data/Step1/expresion_matrix.rds")
sen = readRDS("Processed_Data/Step1/sensitivity_matrix.rds")
dim(GE)
dim(sen)

#loop across drugs--------------------------------------
N_itration = 6

N_drugs = ncol(sen)
mean_mse = rep(0,N_drugs)
sd_mse = rep(0,N_drugs)
mean_corr = rep(0,N_drugs)
sd_corr = rep(0,N_drugs)

for (i in 1:N_drugs){
  print("The drug number is:")
  print(i)
  
  X = GE[!is.na(sen[,i]),]
  dim(X)
  y = sen[!is.na(sen[,i]),i]
  length(y)
  
  Corr = cor(X,y)
  high_corr = order(Corr,decreasing = TRUE)
  X = X[,high_corr[1:200]]
  
  
  # Cross validation loop
  
  clusterExport(cl, c("X","y","i"))
  clusterEvalQ(cl, c(library(caTools),library(randomForest),
                     source("RandomForest_Func.R")))
  
  RepLoop = function(j){
    
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
    mse = mean((ytest-y_pred)^2)
    corr = cor(ytest,y_pred)
    result = data.frame(mse = mse,corr = corr)
    #print(mse)[j]
    #print(corr)[j]
    return(result)
    }
  
  result = parLapply(cl, sapply(1:N_itration, list), RepLoop) 
  
  Result = data.frame()
  for (k in 1:N_itration)
    Result = rbind(Result, result[[k]])

  
  mean_mse[i] = mean(Result[,1])
  sd_mse[i] = sd(Result[,1])
  mean_corr[i] = mean(Result[,2])
  sd_corr[i] = sd(Result[,2])
}
stopCluster(cl)
Results = cbind(mean_mse = mean_mse,sd_mse = sd_mse,
                 mean_corr = mean_corr,sd_corr = sd_corr)
#saveRDS(Results,"Result_All_Drugs.rds")

cor_sort = data.frame(sort(Results[,3],decreasing = TRUE))

