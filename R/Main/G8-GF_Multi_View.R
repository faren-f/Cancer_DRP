rm(list=ls())
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

source("F4-DoRothEA.R")
#source("F9-decoupleR.R")
source("F12-Drug_Targets.R")
library(parallel)

no_cores = detectCores()
cl = makeCluster(no_cores-2)

res_drugs = readRDS("Processed_data/S0/Result_All_Drugs.rds")
order_drugs = data.frame(order = order(res_drugs[,3],decreasing = TRUE))

sen = readRDS("Processed_Data/S1/sensitivity_matrix.rds")
GE = readRDS("Processed_Data/S1/expresion_matrix.rds")

# L1000
l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
O1 = GE[,colnames(GE)%in%l1000_genes]

# Transcription Factors
#O2 = DoRothEA(X = GE)

# Tissue types
O3 = readRDS("Processed_data/S19/sample_tissue_types.rds")

# decoupleR
#O4 = decoupleR(X = GE, method = "gsva")

Results = c()
for (i in order_drugs[3:20,]){
  print(paste0("The drug number is: ", as.character(i)))
  
  # Drug target
  # DTs = Drug_Targets(X= GE)
  # O5 = GE[, DTs[[i]]]
  # 
  # if (length(DTs[[i]])==1){
  #   D = 1
  #   
  # }else if(length(DTs[[i]])<1){
  #   D = 0
  #   
  #   }else{
  #     D = ncol(O5)
  #   }
  
  
  # concatenate all omics data
  omics = GE
  index = rep(1,ncol(GE))
  #index = c(rep(1,ncol(O1)),rep(2,ncol(O2)),rep(3,D))
  
  X = omics[!is.na(sen[,i]),]
  y = sen[!is.na(sen[,i]),i]
  
  Mean_X = apply(X,2,mean)
  STD_X = apply(X,2,sd)
  X = (X-Mean_X)/STD_X
  
  clusterExport(cl, c("X","y","i","index"))
  clusterEvalQ(cl, c(library(caTools),source("F7-RandomForest.R"),
                     library(keras),library(tensorflow),source("F8-MLP.R"),
                     source("F10-Ridge.R"),source("F11-SGL.R")))
  
  RepLoop = function(j){
    
    sample = sample.split(y, SplitRatio = .8)
    
    Xtrain = subset(X, sample == TRUE)
    Xtest  = subset(X, sample == FALSE)
    ytrain = subset(y, sample == TRUE)
    ytest  = subset(y, sample == FALSE)
    
    # Normalization-------------------------------------------------------------
    #Xtrain normalization
    #Mean_X = apply(Xtrain,2,mean)
    #STD_X = apply(Xtrain,2,sd)
    #Xtrain = (Xtrain-Mean_X)/STD_X
    
    # Xtest normalization
    #Xtest = (Xtest-Mean_X)/STD_X
    
    #for when we have tissue types
    # Mean_X = apply(Xtrain[,1:2254],2,mean)
    # STD_X = apply(Xtrain[,1:2254],2,sd)
    # Xtrain_1 = (Xtrain[,1:2254]-Mean_X)/STD_X
    # Xtrain = cbind(Xtrain_1,Xtrain[,2255:ncol(Xtrain)])
    # # Xtest normalization
    # Xtest_1 = (Xtest[,1:2254]-Mean_X)/STD_X
    # Xtest = cbind(Xtest_1,Xtest[2255:ncol(Xtest)])
    
    # Ytrain normalization
    Mean_y = mean(ytrain)
    STD_y = sd(ytrain)
    ytrain_norm = (ytrain-Mean_y)/STD_y
  
    # Models
    y_pred_SGL = My_SGL(ytrain = ytrain_norm ,Xtrain = Xtrain,Xtest = Xtest,index = index)
    y_pred_RF = RandomForest(ytrain = ytrain_norm ,Xtrain = Xtrain,Xtest = Xtest)
    #y_pred_ENet = ENet(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    #y_pred_Lasso = Lasso(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    #y_pred_Ridge = Ridge(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    #y_pred_MLP = MLP(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    
    # y_pred re-normalization
    
    y_pred_SGL = (y_pred_SGL*STD_y)+Mean_y
    y_pred_RF = (y_pred_RF*STD_y)+Mean_y
    #y_pred_ENet = (y_pred_ENet*STD_y)+Mean_y
    #y_pred_Lasso = (y_pred_Lasso*STD_y)+Mean_y
    #y_pred_Ridge = (y_pred_Ridge*STD_y)+Mean_y
    #y_pred_MLP = (y_pred_MLP*STD_y)+Mean_y
    
    # Evaluation
    corr_SGL = cor(ytest,y_pred_SGL)
    corr_RF = cor(ytest,y_pred_RF)
    #corr_ENet = cor(ytest,y_pred_ENet)
    #corr_Lasso = cor(ytest,y_pred_Lasso)
    #corr_Ridge = cor(ytest,y_pred_Ridge)
    #corr_MLP = cor(ytest,y_pred_MLP)
    
    result = data.frame(corr_SGL = corr_SGL, 
                        corr_RF = corr_RF)
                        #corr_ENet = corr_ENet,
                        #corr_Lasso = corr_Lasso,
                        #corr_Ridge = corr_Ridge,
                        #corr_MLP = corr_MLP)
    
    return(result)
  }
  
  N_itration = 12
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


