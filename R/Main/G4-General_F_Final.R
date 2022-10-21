rm(list=ls())
library(parallel)
library(caTools)
library(igraph)
library(dorothea)
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

source("F4-DoRothEA.R")
source("F9-decoupleR.R")

no_cores = detectCores()
cl = makeCluster(no_cores-2)


#Read data--------------------------------------------------
#PRISM
sen = readRDS("Processed_Data/S1/sensitivity_matrix_AUC.rds")
GE = readRDS("Processed_Data/S1/expresion_matrix.rds")
drug_targets = readRDS("Processed_data/S1/drug_targets.rds")
res_drugs = readRDS("Processed_data/Other/Result_All_Drugs.rds")
order_drugs = data.frame(order = order(res_drugs[,3],decreasing = TRUE))

#CCLE
#sen = readRDS("Processed_Data/S15/sensitivity_matrix_Activity_Area.rds")
#GE = readRDS("Processed_Data/S15/expresion_matrix.rds")
drugs = data.frame(colnames(sen))

TF = DoRothEA(X = GE)
TF = t(TF)

Interaction_Network = "STRING"
MyGraph = NA
my_genes = NA
if (Interaction_Network == "STRING"){
  
  GE = readRDS("Processed_data/S9/expresion_matrix_PRISM_STRING.rds")
  STRING_edgelist = readRDS("Processed_data/S9/ppi_STRING_PRISM.rds")
  ppi = rbind(STRING_edgelist$gene_symbol1, STRING_edgelist$gene_symbol2)
  MyGraph = simplify(graph(ppi, directed = FALSE))
  my_genes = V(MyGraph)$name

} else if(Interaction_Network == "Omnipath"){
  
  GE = readRDS("Processed_data/S10/expresion_matrix_PRISM_Omnipath.rds")
  Omnipath_edgelist = readRDS("Processed_data/S10/ppi_Omnipath_PRISM.rds")
  ppi = rbind(Omnipath_edgelist$gene_symbol1, Omnipath_edgelist$gene_symbol2)
  MyGraph = simplify(graph(ppi, directed = FALSE))
  my_genes = V(MyGraph)$name
}
  

#loop across drugs--------------------------------------
N_itration = 6
N_drugs = ncol(sen)
Results = c()
#drug = 1211
#1) i = 1211
i = 1020
print(paste0("The drug number is: ", as.character(i)))
  
X = GE[!is.na(sen[,i]),]
X_TF = TF[!is.na(sen[,i]),]

y = sen[!is.na(sen[,i]),i]

for (N_feat in c(500)){
  print(N_feat)
  # Cross validation loop
  
  clusterExport(cl, c("X","X_TF","y","i","my_genes","MyGraph","Interaction_Network","N_feat"))
  clusterEvalQ(cl, c(source("F1-high_corr.R"), source("F2-mRMR.R"),
                     library(dorothea),source("F4-DoRothEA.R"),
                     library(igraph),source("F3-Infogenes.R"),
                     library(caTools),library(randomForest),source("F7-RandomForest.R"),
                     library(glmnet),library(caret),source ("F6-ENet.R"),
                     library(keras),library(tensorflow),source("F8-MLP.R"),
                     source("F10-Ridge.R")))
  
  RepLoop = function(j){
    
    sample = sample.split(y, SplitRatio = .9)
    
    Xtrain = subset(X, sample == TRUE)
    Xtest  = subset(X, sample == FALSE)
    ytrain = subset(y, sample == TRUE)
    ytest  = subset(y, sample == FALSE)
    
    Xtrain_TF = subset(X_TF, sample == TRUE)
    Xtest_TF  = subset(X_TF, sample == FALSE)
  
    
    # Normalization-------------------------------------------------------------
    # Xtrain normalization
    Mean_X = apply(Xtrain,2,mean)
    STD_X = apply(Xtrain,2,sd)
    Xtrain = (Xtrain-Mean_X)/STD_X
    
    Mean_X_TF = apply(Xtrain_TF,2,mean)
    STD_X_TF = apply(Xtrain_TF,2,sd)
    Xtrain_TF = (Xtrain_TF-Mean_X_TF)/STD_X_TF
    # Xtest normalization
    Xtest = (Xtest-Mean_X)/STD_X
    
    Xtest_TF = (Xtest_TF-Mean_X_TF)/STD_X_TF
    
    # Ytrain normalization
    Mean_y = mean(ytrain)
    STD_y = sd(ytrain)
    ytrain_norm = (ytrain-Mean_y)/STD_y
    
    # Feature selection---------------------------------------------------------
    FS_method_set = c("high_corr", "mRMR","DoRothEA")
    mse_RF = c()
    corr_RF = c()
    mse_Ridge = c()
    corr_Ridge = c()
    for (FS_method in FS_method_set){
      
      if (FS_method == "infogenes"){
        Xtrain_r = Infogenes(Xtrain,ytrain,MyGraph,my_genes,N_feat=N_feat)
        Xtest_r = Xtest[,colnames(Xtrain_r)]
        
      }else if(FS_method == "high_corr"){
        Xtrain_r = high_corr(Xtrain,ytrain,N_feat = N_feat)
        Xtest_r = Xtest[,colnames(Xtrain_r)]
        
      }else if(FS_method == "mRMR"){
        # Xtrain_r = mRMR(Xtrain, ytrain, N_feat = N_feat, alpha=1, do.plot = FALSE)
        # Xtest_r = Xtest[,colnames(Xtrain_r)]
        Xtrain_r = Xtrain
        Xtest_r = Xtest
        
      }else if(FS_method == "DoRothEA"){
        Xtrain_r = high_corr(Xtrain_TF,ytrain,N_feat = N_feat)
        Xtest_r = Xtest_TF[,colnames(Xtrain_r)]
        # Xtrain_r = Xtrain_TF
        # Xtest_r = Xtest_TF
      }
      

      # Models
      
      y_pred_RF = RandomForest(ytrain = ytrain ,Xtrain = Xtrain_r, Xtest = Xtest_r)
      
      y_pred_Ridge = Ridge(ytrain = ytrain ,Xtrain = Xtrain_r, Xtest = Xtest_r)
      
      # y_pred re-normalization
      y_pred_RF = (y_pred_RF*STD_y)+Mean_y
      y_pred_Ridge = (y_pred_Ridge*STD_y)+Mean_y
      
      # Evaluation
      mse_RF = c(mse_RF, mean((ytest-y_pred_RF)^2))
      corr_RF = c(corr_RF, cor(ytest,y_pred_RF))
      
      mse_Ridge = c(mse_Ridge, mean((ytest-y_pred_Ridge)^2))
      corr_Ridge = c(corr_Ridge, cor(ytest,y_pred_Ridge))
    }
    
    # result = data.frame(
    #   mse_RF_HC = mse_RF[1], mse_RF_IG = mse_RF[2], mse_RF_mRMR = mse_RF[3],
    #   cor_RF_HC = corr_RF[1], cor_RF_IG = corr_RF[2], cor_RF_mRMR = corr_RF[3],
    #   mse_Ridge_HC = mse_Ridge[1], mse_Ridge_IG = mse_Ridge[2], mse_Ridge_mRMR = mse_Ridge[3],
    #   cor_Ridge_HC = corr_Ridge[1], cor_Ridge_IG = corr_Ridge[2], cor_Ridge_mRMR = corr_Ridge[3]
    #   )

    result = data.frame(
      cor_RF_HC = corr_RF[1], cor_RF_IG = corr_RF[2], 
      cor_RF_mRMR = corr_RF[3],cor_RF_Do = corr_RF[4],
      cor_Ridge_HC = corr_Ridge[1], cor_Ridge_IG = corr_Ridge[2], 
      cor_Ridge_mRMR = corr_Ridge[3],cor_Ridge_Do = corr_Ridge[4]
    )   
    
    return(result)
  }
  
  result = parLapply(cl, sapply(1:N_itration, list), RepLoop) 
  
  Result = data.frame()
  for (k in 1:N_itration){
    Result = rbind(Result, result[[k]])
  }
  
  Results = rbind(Results, apply(Result, 2, mean))
  print(apply(Result, 2, mean))
}
Results = data.frame(Results)
stopCluster(cl)


# plot(Results$cor_RF_HC, type = "l", ylim = c(.4,.65))
# par(new = TRUE)
# plot(Results$cor_RF_IG, type = "l", ylim = c(.4,.65), col = "red")
# par(new = TRUE)
# plot(Results$cor_RF_mRMR, type = "l", ylim = c(.4,.65), col = "blue")

