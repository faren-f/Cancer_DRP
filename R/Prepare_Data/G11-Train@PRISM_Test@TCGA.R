rm(list=ls())

source("F7-RandomForest.R")
source("F6-ENet.R")
source("F8-MLP.R")
source("F10-Ridge.R")
source("F11-SGL.R")
source("F13-Lasso.R")

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
sen_PRISM = readRDS("Processed_data/S23/sensitivity_matrix_PRISM_with@TCGA@drugs.rds")
res_TCGA = readRDS("Processed_data/S23/Drug_response_matrix_TCGA.rds")
GE = readRDS("Processed_data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")
GE_TCGA = readRDS("Processed_data/S23/expresion_matrix_TCGA.rds")
GE_TCGA[is.na(GE_TCGA)] = 0
r = apply(GE_TCGA,2,quantile,prob=0.75)
hist(r)
GE_TCGA = GE_TCGA[,-which(r==0)]
GE = GE[,-which(r==0)]

N_drug = ncol(sen_PRISM)
Results = c()
for (i in 1){
  i = 1
  print(paste0("The drug number is: ", as.character(i)))
  
  Xtrain = GE[!is.na(sen_PRISM[,i]),]
  ytrain = sen_PRISM[!is.na(sen_PRISM[,i]),i]

  Mean_Xtrain = apply(Xtrain,2,mean)
  STD_Xtrain = apply(Xtrain,2,sd)
  Xtrain = (Xtrain-Mean_Xtrain)/STD_Xtrain
  hist(Xtrain)
  
  Xtest = GE_TCGA[!is.na(res_TCGA[,i]),]
  ytest = res_TCGA[!is.na(res_TCGA[,i]),i]
  
  Mean_Xtest = apply(Xtest,2,mean)
  STD_Xtest = apply(Xtest,2,sd)
  Xtest = (Xtest-Mean_Xtest)/STD_Xtest
  hist(Xtest)
  
  source("F15-Feature_Selection_PRISM@TCGA.R")
  selected_features = c("Landmark_genes")
  Omics_List = Feature_Selection(selected_features,GE = Xtrain ,GE_test = Xtest)
  Xtrain = Omics_List[[1]]
  index = Omics_List[[2]]
  Xtest = Omics_List[[3]]
  
  # Ytrain normalization
  Mean_ytrain = mean(ytrain)
  STD_ytrain = sd(ytrain)
  ytrain = (ytrain-Mean_ytrain)/STD_ytrain
  
  # Models
  y_pred_SGL = My_SGL(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest,index = index)
  #y_pred_RF = RandomForest(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
  #y_pred_ENet = ElasticNet(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
  #y_pred_Lasso = Lasso(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
  y_pred_Ridge = Ridge(ytrain = ytrain ,Xtrain = Xtrain, Xtest = Xtest)
  #y_pred_MLP = MLP(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
  
  # Evaluation
  #corr_SGL = cor(ytest,y_pred_SGL)
  #corr_RF = cor(ytest,y_pred_RF)
  #corr_ENet = cor(ytest,y_pred_ENet)
  #corr_Lasso = cor(ytest,y_pred_Lasso)
  corr_Ridge = cor(ytest,y_pred_Ridge)
  t.test(y_pred_SGL[ytest<=2], y_pred_SGL[ytest>=3])
  plot(ytest,y_pred_Ridge)
  plot(ytest,y_pred_SGL)
  #corr_MLP = cor(ytest,y_pred_MLP)
  
  result = data.frame(corr_Ridge = corr_Ridge)
                      #corr_SGL = corr_SGL,
                      #corr_RF = corr_RF,
                      #corr_ENet = corr_ENet,
                      #corr_Lasso = corr_Lasso,
                      #corr_MLP = corr_MLP)
  
  
  Result_mean = apply(Result, 2, mean)
  Result_sd = apply(Result, 2, sd)
  print(Result_mean)
  
  Results = rbind(Results, c(Result_mean, Result_sd))
  
}



