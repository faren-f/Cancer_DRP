rm(list=ls())

library(caTools)
source("F14-Feature_Selection.R")
source("F10-Ridge.R")

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
GE_PRISM = readRDS("Processed_Data/S1/expresion_matrix.rds")
sen_PRISM = readRDS("Processed_data/S1/sensitivity_matrix_AUC.rds")

GE_GDSC = readRDS("Processed_Data/S41/GDSC_expresion_matrix.rds")
sen_GDSC = readRDS("Processed_Data/S41/GDSC_sensitivity_matrix_AUC.rds")

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

Mean_Corr = c()
STD_Corr = c()

Mean_MSE = c()
STD_MSE = c()

for (i in 1:ncol(sen)){             # drug loop
  print(i)
  Corr = c()
  MSE = c()
  
  for (j in 1:1){           # repeat loop
    print(j)
    Xtrain = GE_PRISM[!is.na(sen_PRISM[,i]),]
    ytrain = sen_PRISM[!is.na(sen_PRISM[,i]),i]
    
    Xtest = GE_GDSC[!is.na(sen_GDSC[,i]),]
    ytest = sen_GDSC[!is.na(sen_GDSC[,i]),i]
    
    # X normalization
    Xtrain = scale(Xtrain)
    Xtest = scale(Xtest)
    
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
  
  Mean_Corr = rbind(Mean_Corr, mean(Corr))
  STD_Corr = rbind(STD_Corr, sd(Corr))
  
  Mean_MSE = rbind(Mean_MSE, mean(MSE))
  STD_MSE = rbind(STD_MSE, sd(MSE))
  
}

rownames(Mean_Corr) = colnames(sen)

rownames(STD_Corr) = colnames(sen)

rownames(Mean_MSE) = colnames(sen)

rownames(STD_MSE) = colnames(sen)

print(Mean_Corr)

