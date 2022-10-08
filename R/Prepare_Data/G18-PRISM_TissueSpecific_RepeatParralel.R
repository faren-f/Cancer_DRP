rm(list=ls())

library(parallel)
no_cores = detectCores()
cl = makeCluster(no_cores-2)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
source("F14-Feature_Selection.R")

Cellline_Tissue = readRDS("Processed_data/S29/Cellline_Tissue.rds")
Tissue_Clusters = readRDS("Processed_data/S29/Tissue_Clusters.rds")
table(Tissue_Clusters)

sen = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")
GE = readRDS("Processed_Data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")

# Seperating Tissues whithin the first cluster
#Tissues_k1 = which(Cellline_Tissue[,1]==7|Cellline_Tissue[,1]==21) # methotrexate
#Tissues_k1 = which(Cellline_Tissue[,1]==11|Cellline_Tissue[,1]==13|Cellline_Tissue[,1]==15) #dacarbazine

#Cellline_Tissue_k1 = Cellline_Tissue[Tissues_k1,]

#sen = sen[rownames(sen) %in% rownames(Cellline_Tissue_k1),]
#GE = GE[rownames(GE) %in% rownames(Cellline_Tissue_k1),]

# Tissues_k1 = which(Tissue_Clusters==2)
# Cellline_Tissue_k1 = Cellline_Tissue[Cellline_Tissue[,1] %in% Tissues_k1,]
# 
# sen = sen[rownames(Cellline_Tissue_k1),]
# GE = GE[rownames(Cellline_Tissue_k1),]

selected_features = c("Landmark_genes")
Omics_List = Feature_Selection(selected_features,GE)
omics = Omics_List[[1]]
index = Omics_List[[2]]

#omics = GE
#index = rep(1,ncol(GE))

N_drug = ncol(sen)
Results = c()
i=12
for (i in 12){
  print(paste0("The drug number is: ", as.character(i)))
  
  X = omics[!is.na(sen[,i]),]
  y = sen[!is.na(sen[,i]),i]
  
  #X = scale(X)
  
  # Ytrain normalization
  #Mean_y = mean(y)
  #STD_y = sd(y)
  #y = (y-Mean_y)/STD_y
  
  source("F10-Ridge.R")
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
    #y_pred_MLP = MLP(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    
    
    # Evaluation
    #corr_SGL = cor(ytest,y_pred_SGL)
    #corr_RF = cor(ytest,y_pred_RF)
    #corr_ENet = cor(ytest,y_pred_ENet)
    #corr_Lasso = cor(ytest,y_pred_Lasso)
    corr_Ridge = cor(ytest,y_pred_Ridge,method = "pearson")
    #corr_Ridge
    #corr_MLP = cor(ytest,y_pred_MLP)
    #plot(ytest,y_pred_Ridge,xlim = c(0,1.4), ylim = c(0,1.4))
    result = data.frame(corr_Ridge = corr_Ridge)
    #corr_SGL = corr_SGL,
    #corr_RF = corr_RF,
    #corr_ENet = corr_ENet,
    #corr_Lasso = corr_Lasso,
    #corr_MLP = corr_MLP)
    
    return(result)
  }
  
  N_itration = 50
  result = parLapply(cl, sapply(1:N_itration, list), RepLoop) 
  
  Result = data.frame()
  for (k in 1:N_itration){
    Result = rbind(Result, result[[k]])
  }
  
  
  Result_mean = apply(Result, 2, mean)
  Result_sd = apply(Result, 2, sd)
  print(Result_mean)
  print(Result_sd)
  
  
  Results = rbind(Results, c(Result_mean, Result_sd))
  
}
stopCluster(cl)




