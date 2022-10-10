rm(list=ls())

library(parallel)
no_cores = detectCores()
cl = makeCluster(no_cores-2)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
sen = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")
GE = readRDS("Processed_Data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")


source("F14-Feature_Selection.R")
selected_features = c("Landmark_genes")
Omics_List = Feature_Selection(selected_features,GE)
omics = Omics_List[[1]]
index = Omics_List[[2]]

N_drug = ncol(sen)
Results = c()
for (i in 5){
  print(paste0("The drug number is: ", as.character(i)))
  
  X = omics[!is.na(sen[,i]),]
  y = sen[!is.na(sen[,i]),i]
  #boxplot(y)
  
  # U = 1.5*IQR(y)+median(y)
  # L = median(y)-1.5*IQR(y)
  # X = X[y<U & y>L,]
  # y = y[y<U & y>L]
  
  
  X = t(scale(t(X)))
  
  
  clusterExport(cl, c("X","y","i","index"))
  clusterEvalQ(cl, c(library(caTools),source("F7-RandomForest.R"),
                     source("F6-ENet.R"),source("F8-MLP.R"),source("F10-Ridge.R"),
                     source("F11-SGL.R"),source("F13-Lasso.R")))
  
  RepLoop = function(j){
    
    sample = sample.split(y, SplitRatio = .9)
    
    Xtrain = subset(X, sample == TRUE)
    Xtest  = subset(X, sample == FALSE)
    ytrain = subset(y, sample == TRUE)
    ytest  = subset(y, sample == FALSE)
    
    # a = hist(ytrain)
    # InvCount = 1/a[["counts"]]
    # 
    # W = rep(0,length(ytrain))
    # b = a$breaks
    # for(l in 1:length(ytrain)){
    #   for(k in 2:length(b)){
    #     if((b[k-1] < ytrain[l]) & (ytrain[l] < b[k])){
    #     W[l] = InvCount[k-1]
    #     }
    #   }
    # }
    # hist(W)
    
    #W = rep(1,length(ytrain))
    
    # Models
    #y_pred_SGL = My_SGL(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest,index = index)
    #y_pred_RF = RandomForest(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    #y_pred_ENet = ElasticNet(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    #y_pred_Lasso = Lasso(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    #y_pred_Ridge = Ridge(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest, weight = W)
    #y_pred_MLP = MLP(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)

    
    
    # Evaluation
    #corr_SGL = cor(ytest,y_pred_SGL)
    #corr_RF = cor(ytest,y_pred_RF)
    #corr_ENet = cor(ytest,y_pred_ENet)
    #corr_Lasso = cor(ytest,y_pred_Lasso)
    corr_Ridge = cor(ytest,y_pred_Ridge,method = "pearson")
    #corr_MLP = cor(ytest,y_pred_MLP)
    
    
    result = data.frame(corr_Ridge = corr_Ridge)
    #corr_SGL = corr_SGL,
    #corr_RF = corr_RF,
    #corr_ENet = corr_ENet,
    #corr_Lasso = corr_Lasso,
    #corr_MLP = corr_MLP)
    
    return(result)
  }
  
  N_itration = 1
  result = parLapply(cl, sapply(1:N_itration, list), RepLoop) 
  
  Result = data.frame()
  for (k in 1:N_itration){
    Result = rbind(Result, result[[k]])
  }
  
  Result_mean = apply(Result, 2, mean)
  #Result_sd = apply(Result, 2, sd)
  print(Result_mean)
  #print(Result_sd)
  
  
  Results = rbind(Results, Result_mean)
  
}
stopCluster(cl)



