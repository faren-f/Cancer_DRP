rm(list=ls())

library(parallel)
no_cores = detectCores()
cl = makeCluster(no_cores-2)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
sen = readRDS("All_Results/sen_PRISM_good_drugs.rds")
TCGA_PRISM_drugs = readRDS("Processed_data/S21/TCGA_PRISM_drugs.rds")
which(colnames(sen) %in% TCGA_PRISM_drugs)


source("F14-Feature_Selection.R")
selected_features = c("TF_DoRothEA")
Omics_List = Feature_Selection(selected_features)
omics = Omics_List[[1]]
index = Omics_List[[2]]

N_drug = ncol(sen)
Results = c()
#4,5,14,16,17,31,39,48,75,78
for (i in 14){
  print(paste0("The drug number is: ", as.character(i)))

  X = omics[!is.na(sen[,i]),]
  y = sen[!is.na(sen[,i]),i]
  
  Mean_X = apply(X,2,mean)
  STD_X = apply(X,2,sd)
  X = (X-Mean_X)/STD_X
  
  # Ytrain normalization
  Mean_y = mean(y)
  STD_y = sd(y)
  y = (y-Mean_y)/STD_y
  

  clusterExport(cl, c("X","y","i","index"))
  clusterEvalQ(cl, c(library(caTools),source("F7-RandomForest.R"),
                     source("F6-ENet.R"),source("F8-MLP.R"),source("F10-Ridge.R"),
                     source("F11-SGL.R"),source("F13-Lasso.R")))
  
  RepLoop = function(j){
    
    sample = sample.split(y, SplitRatio = .8)
    
    Xtrain = subset(X, sample == TRUE)
    Xtest  = subset(X, sample == FALSE)
    ytrain = subset(y, sample == TRUE)
    ytest  = subset(y, sample == FALSE)
    
    
    # Models
    #y_pred_SGL = My_SGL(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest,index = index)
    #y_pred_RF = RandomForest(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    #y_pred_ENet = ElasticNet(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    #y_pred_Lasso = Lasso(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    y_pred_Ridge = Ridge(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    y_pred_MLP = MLP(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    
    # Evaluation
    #corr_SGL = cor(ytest,y_pred_SGL)
    #corr_RF = cor(ytest,y_pred_RF)
    #corr_ENet = cor(ytest,y_pred_ENet)
    #corr_Lasso = cor(ytest,y_pred_Lasso)
    corr_Ridge = cor(ytest,y_pred_Ridge)
    corr_MLP = cor(ytest,y_pred_MLP)
    
    result = data.frame(corr_Ridge = corr_Ridge,
                        #corr_SGL = corr_SGL,
                        #corr_RF = corr_RF,
                        #corr_ENet = corr_ENet,
                        #corr_Lasso = corr_Lasso,
                        corr_MLP = corr_MLP)
    
    return(result)
  }
  
  N_itration = 6
  result = parLapply(cl, sapply(1:N_itration, list), RepLoop) 
  
  Result = data.frame()
  for (k in 1:N_itration){
    Result = rbind(Result, result[[k]])
  }
  
  
  Result_mean = apply(Result, 2, mean)
  Result_sd = apply(Result, 2, sd)
  print(Result_mean)
  
  Results = rbind(Results, c(Result_mean, Result_sd))
  
}
stopCluster(cl)
#c = colnames(sen)
#rownames(Results) = c[order_drugs[3:20,]]
#saveRDS(Results,"All_Results/.rds")
#all = readRDS("All_Results/SGL_RF@L1000_TF@12run.rds")
#a = data.frame(colnames(sen))



