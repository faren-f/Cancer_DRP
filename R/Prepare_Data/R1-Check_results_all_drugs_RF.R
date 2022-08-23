rm(list=ls())
library(parallel)
library(caTools)

no_cores = detectCores()
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

#Read data--------------------------------------------------
#PRISM
sen = readRDS("Processed_Data/S1/sensitivity_matrix_AUC.rds")
GE = readRDS("Processed_Data/S9/expresion_matrix_PRISM_STRING.rds")


#loop across drugs--------------------------------------
N_itration = 10
N_drugs = ncol(sen)
Results = c()

for (i in 1:N_drugs){
  print(paste0("The drug number is: ", as.character(i)))

  X = GE[!is.na(sen[,i]),]
  y = sen[!is.na(sen[,i]),i]
  
  # Normalization-------------------------------------------------------------
  # Xtrain normalization
  Mean_X = apply(X,2,mean)
  STD_X = apply(X,2,sd)
  X = (X-Mean_X)/STD_X
  
  # Ytrain normalization
  Mean_y = mean(y)
  STD_y = sd(y)
  y = (y-Mean_y)/STD_y
  
  # Cross validation loop
  cl = makeCluster(no_cores-1)
  clusterExport(cl, c("X","y","i"))
  clusterEvalQ(cl, c(library(caTools),source("F7-RandomForest.R")))
  
  RepLoop = function(j){
    
    sample = sample.split(y, SplitRatio = .9)
    
    Xtrain = subset(X, sample == TRUE)
    Xtest  = subset(X, sample == FALSE)
    ytrain = subset(y, sample == TRUE)
    ytest  = subset(y, sample == FALSE)
    
    # Models
    y_pred_RF = RandomForest(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    
    # Evaluation
    mse_RF = mean((ytest-y_pred_RF)^2)
    corr_RF = cor(ytest,y_pred_RF)
    result = data.frame(corr_RF = corr_RF)
    
    return(result)
  }

  result = parLapply(cl, sapply(1:N_itration, list), RepLoop) 
  
  Result = data.frame()
  for (k in 1:N_itration){
    Result = rbind(Result, result[[k]])
  }
  
  r = c(mean(Result$corr_RF),sd(Result$corr_RF))
  
  print(r[1])
  Results = rbind(Results, r)

  stopCluster(cl)
}
rownames(Results) = colnames(sen)
saveRDS(Results,"Processed_Data/Result_RF_All_Drugs.rds")

