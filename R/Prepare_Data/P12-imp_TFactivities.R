rm(list=ls())

library(ROCR)
source("F18-Combat_Normalization.R")
library(glmnet)
library(caret)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

dR_PRISM = read.table("Processed_data/S33/gsea2_PRISM.csv",sep = ",",header = TRUE, row.names = 1)
dR_TCGA = read.table("Processed_data/S33/gsea2_TCGA.csv",sep = ",",header = TRUE, row.names = 1)


sen_PRISM = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")
res_TCGA = readRDS("Processed_data/Other/Res_TCGA_24_Drugs.rds")

# Remove genes whose Q3 is zero
q3_genes = apply(dR_TCGA,2,quantile,prob=0.75)
if(sum(q3_genes==0)>0){
  dR_TCGA = dR_TCGA[,-which(q3_genes==0)]
  dR_PRISM = dR_PRISM[,-which(q3_genes==0)]
}

SigDrugs = c("etoposide","paclitaxel","leucovorin", 
             "ifosfamide", "gemcitabine",
             "cisplatin", "vinblastine")

result = c()
Wilcox_Test_AllDrugs = c()
T_Test_AllDrugs = c()

for (i in SigDrugs){
  print(paste0("The drug number is: ", i))
  
  Xtrain = dR_PRISM[!is.na(sen_PRISM[,i]),]
  ytrain = sen_PRISM[!is.na(sen_PRISM[,i]),i]
  
  Xtest = dR_TCGA[!is.na(res_TCGA[,i]),]
  ytest = res_TCGA[!is.na(res_TCGA[,i]),i]
  
  X_Normalization = Combat_Scale(Xtrain,Xtest)
  
  Xtrain = X_Normalization[[1]]
  Xtest = X_Normalization[[2]]
  
  # Ytrain normalization
  ytrain = scale(ytrain)
  ytrain = ytrain[,1]
  
  train_data = cbind(Xtrain,ytrain)
  control = trainControl(method = "repeatedcv",
                         number = 5,
                         repeats = 5,
                         verboseIter = FALSE)
  
  tune = expand.grid(alpha = 0,lambda = seq(0.01,5,by = 0.01))
  
  model = caret::train(ytrain ~., data = train_data,
                       method = "glmnet",
                       metric="RMSE",
                       allowParallel = TRUE,
                       tuneGrid = tune,
                       trControl = control)
  y_pred = predict(model,Xtest)
  Beta = as.matrix(coef(model$finalModel, model$bestTune$lambda))

  Beta_Null = c()
  Wilcox_Test = c()
  T_test = c()
  for(j in 1:100){
  
    ytrain_perm = sample(ytrain)
    
    train_data = cbind(Xtrain,ytrain_perm)
    control = trainControl(method = "repeatedcv",
                           number = 5,
                           repeats = 5,
                           verboseIter = FALSE)
    
    tune = expand.grid(alpha = 0,lambda = seq(0.01,5,by = 0.05))
    
    model = caret::train(ytrain_perm ~., data = train_data,
                         method = "glmnet",
                         metric="RMSE",
                         allowParallel = TRUE,
                         tuneGrid = tune,
                         trControl = control)
    y_pred = predict(model,Xtest)
    beta = as.matrix(coef(model$finalModel, model$bestTune$lambda))
    Beta_Null = cbind(Beta_Null, beta)
  }
  for(k in 1:nrow(Beta_Null)){
    W = wilcox.test(abs(Beta_Null[k,]), Beta[k], alternative = "greater")$p.value
    Wilcox_Test = c(Wilcox_Test, W)
    
    t_test = t.test(abs(Beta_Null[k,]), mu = Beta[k], alternative = "greater")$p.value
    T_test = c(T_test, t_test)
  }
  Wilcox_Test_AllDrugs = cbind(Wilcox_Test_AllDrugs,Wilcox_Test)
  T_Test_AllDrugs = cbind(T_Test_AllDrugs, T_test)
  
}

#sum(Wilcox_Test_AllDrugs[,1]<0.05)
#which(Wilcox_Test_AllDrugs[,1]<0.05)


