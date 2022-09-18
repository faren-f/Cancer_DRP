rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

source("F7-RandomForest.R")
source("F6-ENet.R")
source("F8-MLP.R")
source("F10-Ridge.R")
source("F11-SGL.R")
source("F13-Lasso.R")
source("F18-Combat_Normalization.R")


sen_PRISM = readRDS("Processed_data/S23/sensitivity_matrix_PRISM_with@TCGA@drugs.rds")
res_TCGA = readRDS("Processed_data/S24/Drug_response_TCGA_binarized.rds")

GE = readRDS("Processed_data/S25/expresion_matrix_Combat_PRISM_with@TCGA@genes.rds")
GE_TCGA = readRDS("Processed_data/S25/expresion_matrix_Combat_TCGA.rds")

N_drug = ncol(sen_PRISM)
drugs = data.frame(colnames(sen_PRISM))
Results = c()

for (i in 1:N_drug){
  
  print(paste0("The drug number is: ", as.character(i)))
  
  Xtrain = GE[!is.na(sen_PRISM[,i]),]
  ytrain = sen_PRISM[!is.na(sen_PRISM[,i]),i]
  
  Xtest = GE_TCGA[!is.na(res_TCGA[,i]),]
  ytest = res_TCGA[!is.na(res_TCGA[,i]),i]
  
  length(ytest)
  if(length(ytest)>9){
    
    source("F15-Feature_Selection_PRISM@TCGA.R")
    selected_features = "Landmark_genes"
    Omics_List = Feature_Selection(selected_features,GE = Xtrain ,GE_test = Xtest)
    Xtrain = Omics_List[[1]]
    index = Omics_List[[2]]
    Xtest = Omics_List[[3]]
    
    # Ytrain normalization
    #Mean_ytrain = mean(ytrain)
    #STD_ytrain = sd(ytrain)
    #ytrain = (ytrain-Mean_ytrain)/STD_ytrain
    
    # Models
    #y_pred_SGL = My_SGL(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest,index = index)
    #y_pred_RF = RandomForest(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    #y_pred_ENet = ElasticNet(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    #y_pred_Lasso = Lasso(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    y_pred_Ridge = Ridge(ytrain = ytrain ,Xtrain = Xtrain, Xtest = Xtest)
    #y_pred_Ridge = MLP(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    
    # Evaluation
    #corr_SGL = cor(ytest,y_pred_SGL)
    #corr_RF = cor(ytest,y_pred_RF)
    #corr_ENet = cor(ytest,y_pred_ENet)
    #corr_Lasso = cor(ytest,y_pred_Lasso)
    corr_Ridge = cor(ytest , y_pred_Ridge)
    
    ttest = t.test(y_pred_Ridge[ytest==1], y_pred_Ridge[ytest==2], alternative="greater")$p.value
    Ranksum = wilcox.test(y_pred_Ridge[ytest==1], y_pred_Ridge[ytest==2], alternative ="greater")$p.value
    
  } else{
    corr_Ridge = 0
    ttest = 1
    Ranksum = 1
  }
  plot(ytest,y_pred_Ridge)
  boxplot(y_pred_Ridge[ytest==1], y_pred_Ridge[ytest==2])
  #plot(ytest,y_pred_SGL)
  #corr_MLP = cor(ytest,y_pred_MLP)
  result = data.frame(corr_Ridge = corr_Ridge, ttest=ttest, Ranksum = Ranksum)
  #corr_SGL = corr_SGL,
  #corr_RF = corr_RF,
  #corr_ENet = corr_ENet,
  #corr_Lasso = corr_Lasso,
  #corr_MLP = corr_MLP)
  Results = rbind(Results, result)
  print(Results)
}

sum(Results$Ranksum<0.05)
which(Results$Ranksum<0.05)

