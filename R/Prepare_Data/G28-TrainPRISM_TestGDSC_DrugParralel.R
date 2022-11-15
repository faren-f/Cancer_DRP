rm(list=ls())

library(parallel)
no_cores = detectCores()
cl = makeCluster(no_cores-2)

library(caTools)
source("F14-Feature_Selection.R")
source("F10-Ridge.R")

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
GE_PRISM = readRDS("Processed_Data/S1/expresion_matrix.rds")
sen_PRISM = readRDS("Processed_data/S1/sensitivity_matrix_AUC.rds")

GE_GDSC = readRDS("Processed_Data/S39/GDSC_expresion_matrix.rds")
sen_GDSC = readRDS("Processed_Data/S39/GDSC_sensitivity_matrix_AUC.rds")

intersected_genes = intersect(colnames(GE_PRISM),colnames(GE_GDSC))
GE_PRISM = GE_PRISM[,intersected_genes]
GE_GDSC = GE_GDSC[,intersected_genes]

intersect_drugs = intersect(colnames(sen_PRISM), colnames(sen_GDSC))
sen_PRISM = sen_PRISM[,intersect_drugs]
sen_GDSC = sen_GDSC[,intersect_drugs]

selected_features = c("Landmark_genes")
Omics_List = Feature_Selection(selected_features,GE_PRISM)
GE_PRISM = Omics_List[[1]]
GE_GDSC = GE_GDSC[,colnames(GE_PRISM)]


clusterExport(cl, c("GE_PRISM","GE_GDSC","sen_PRISM","sen_GDSC"))
clusterEvalQ(cl, c(source("F10-Ridge.R"),source("F18-Combat_Normalization.R")))

DrugLoop = function(i){
  
  Corr = c()
  MSE = c()
  
  for (j in 1:1){           # repeat loop
    print(j)
    Xtrain = GE_PRISM[!is.na(sen_PRISM[,i]),]
    ytrain = sen_PRISM[!is.na(sen_PRISM[,i]),i]
    
    Xtest = GE_GDSC[!is.na(sen_GDSC[,i]),]
    ytest = sen_GDSC[!is.na(sen_GDSC[,i]),i]
    
    
    X_Normalization = Combat_Scale(Xtrain,Xtest)
    Xtrain = X_Normalization[[1]]
    Xtest = X_Normalization[[2]]
    
    # Ytrain normalization
    ytrain = scale(ytrain)
    ytrain = ytrain[,1]
    
    # Models
    y_pred = Ridge(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    
    # Evaluation
    corr = cor(ytest,y_pred)
    Corr = c(Corr, corr)
    
    mse =  mean((ytest - y_pred)^2)
    MSE = c(MSE, mse)
  }
  
  result = data.frame(Mean_Corr = mean(Corr), STD_Corr = sd(Corr),
                     Mean_MSE = mean(MSE), STD_MSE = sd(MSE))
  return(result)
}

N_drug = ncol(sen_PRISM)
result = parLapply(cl, sapply(1:N_drug, list), DrugLoop) 

Result = data.frame()
for (k in 1:N_drug){
  Result = rbind(Result, result[[k]])
}

stopCluster(cl)
rownames(Result) = colnames(sen_PRISM)
mean(Result[,1])
sd(Result[,3])
hist(Result[,1])
