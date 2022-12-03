rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

library(ROCR)
source("F18-Combat_Normalization.R")
library(glmnet)
library(caret)

sen_PRISM = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")
res_TCGA = readRDS("Processed_data/Other/Res_TCGA_24_Drugs.rds")

GE_PRISM = readRDS("Processed_data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")
GE_TCGA = readRDS("Processed_data/S23/expresion_matrix_TCGA.rds")

# Remove genes whose Q3 is zero
q3_genes = apply(GE_TCGA,2,quantile,prob=0.75)
if(sum(q3_genes==0)>0){
  GE_TCGA = GE_TCGA[,-which(q3_genes==0)]
  GE_PRISM = GE_PRISM[,-which(q3_genes==0)]
}

SigDrugs = c("vinorelbine","gemcitabine","cisplatin")
result = c()
Wilcox_Test_AllDrugs = c()
T_Test_AllDrugs = c()

for (i in SigDrugs){
  print(paste0("The drug number is: ", i))
  
  Xtrain = GE_PRISM[!is.na(sen_PRISM[,i]),]
  ytrain = sen_PRISM[!is.na(sen_PRISM[,i]),i]
  
  Xtest = GE_TCGA[!is.na(res_TCGA[,i]),]
  ytest = res_TCGA[!is.na(res_TCGA[,i]),i]
  
  X_Normalization = Combat_Scale(Xtrain,Xtest)
  
  Xtrain = X_Normalization[[1]]
  Xtest = X_Normalization[[2]]
  
  source("F15-Feature_Selection_PRISM@TCGA.R")
  selected_features = c("progeny")
  Omics_List = Feature_Selection_PRISM_TCGA(selected_features, Xtrain=Xtrain, Xtest=Xtest)
  Xtrain = Omics_List[[1]]
  index = Omics_List[[2]]
  Xtest = Omics_List[[3]]
  
  # Ytrain normalization
  ytrain = scale(ytrain)
  ytrain = ytrain[,1]
  
  train_data = cbind(Xtrain,ytrain)
  control = trainControl(method = "repeatedcv",
                         number = 5,
                         repeats = 5,
                         verboseIter = FALSE)

  tune = expand.grid(alpha = 0,lambda = seq(0.01,5,by = 0.05))
  
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
    W = wilcox.test(Beta_Null[k,],Beta[k], alternative = "two.sided")$p.value
    Wilcox_Test = c(Wilcox_Test, W)
    
    t_test = t.test(Beta_Null[k,], mu = Beta[k], alternative = "two.sided")$p.value
    T_test = c(T_test, t_test)
  }
  Wilcox_Test_AllDrugs = cbind(Wilcox_Test_AllDrugs,Wilcox_Test)
  T_Test_AllDrugs = cbind(T_Test_AllDrugs, T_test)
  
}
sum(Wilcox_Test_AllDrugs[,1]<0.05)
#which(Wilcox_Test_AllDrugs[,1]<0.05)






