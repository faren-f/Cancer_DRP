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

data_GE = c("GE_PRISM", "GE_GDSC")
data_sen = c("sen_PRISM","sen_GDSC")

Mean_Corr_datasets = c()
STD_Corr_datasets = c()

Mean_MSE_datasets = c()
STD_MSE_datasets = c()

for (d in 1:length(data_GE)){                              # model loop
  GE = get(data_GE[d])
  sen = get(data_sen[d])
  
  Mean_Corr = c()
  STD_Corr = c()
  
  Mean_MSE = c()
  STD_MSE = c()
  
  for (i in 1:ncol(sen)){             # drug loop
    print(i)
    Corr = c()
    MSE = c()
    
    for (j in 1:50){           # repeat loop
      print(j)
      X = GE[!is.na(sen[,i]),]
      y = sen[!is.na(sen[,i]),i]
      
      sample = sample.split(y, SplitRatio = .8)
      
      Xtrain = subset(X, sample == TRUE)
      Xtest  = subset(X, sample == FALSE)
      ytrain = subset(y, sample == TRUE)
      ytest  = subset(y, sample == FALSE)
      
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
  Mean_Corr_datasets = cbind(Mean_Corr_datasets, Mean_Corr)
  STD_Corr_datasets = cbind(STD_Corr_datasets, STD_Corr)
  
  Mean_MSE_datasets = cbind(Mean_MSE_datasets, Mean_MSE)
  STD_MSE_datasets = cbind(STD_MSE_datasets, STD_MSE)
}

rownames(Mean_Corr_datasets) = colnames(sen)
colnames(Mean_Corr_datasets) = c("Mean_Corr_PRISM", "Mean_Corr_GDSC")

rownames(STD_Corr_datasets) = colnames(sen)
colnames(STD_Corr_datasets) = c("STD_Corr_PRISM", "STD_Corr_GDSC")

rownames(Mean_MSE_datasets) = colnames(sen)
colnames(Mean_MSE_datasets) = c("Mean_MSE_PRISM", "Mean_MSE_GDSC")

rownames(STD_MSE_datasets) = colnames(sen)
colnames(STD_MSE_datasets) = c("STD_MSE_PRISM", "STD_MSE_GDSC")


Result = cbind(Mean_Corr_datasets,STD_Corr_datasets,Mean_MSE_datasets,STD_MSE_datasets)
plot(Result[,1],Result[,2], xlim = c(0,0.6), ylim = c(0,0.6))

#saveRDS(Result,"Final_Result/Comparison_PRISM&GDSC/Result_PRISM&GDSC_Comparison.rds")
