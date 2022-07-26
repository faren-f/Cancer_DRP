rm(list=ls())

library(caTools)
library(parallel)
library(igraph)

no_cores = detectCores()
cl = makeCluster(no_cores-2)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

#Read data--------------------------------------------------
sen = readRDS("Processed_Data/Step1/sensitivity_matrix.rds")

Interaction_Network = "STRING"

if (Interaction_Network == "None"){
  GE = readRDS("Processed_Data/Step1/expresion_matrix.rds")
  
} else if(Interaction_Network == "STRING"){
  
  GE = readRDS("Processed_data/Step9/expresion_matrix_PRISM_STRING.rds")
  STRING_edgelist = readRDS("Processed_data/Step9/ppi_STRING_PRISM.rds")
  ppi = rbind(STRING_edgelist$gene_symbol1, STRING_edgelist$gene_symbol2)
  MyGraph = simplify(graph(ppi, directed = FALSE))
  my_genes = V(MyGraph)$name
  
} else if(Interaction_Network == "Ompnipath"){
  
  GE = readRDS("Processed_data/Step10/expresion_matrix_PRISM_Omnipath.rds")
  Omnipath_edgelist = readRDS("Processed_data/Step10/ppi_Omnipath_PRISM.rds")
  ppi = rbind(Omnipath_edgelist$gene_symbol1, Omnipath_edgelist$gene_symbol2)
  MyGraph = simplify(graph(ppi, directed = FALSE))
  my_genes = V(MyGraph)$name
  
} else if(Interaction_Network == "OP_DoRothEA"){
  TF = readRDS("Processed_data/Step13/TF_wmean.rds")
}


omics = "GE"
if (omics =="GE"){
  omic = GE
  
}else if(omics =="TF"){
  omic = TF
  
}else if(omics =="Mu"){
  omic = Mu
  
}else if(omics =="GE+Mu+CNV"){
  omic = c(GE,Mu,CNV)
}


#loop across drugs--------------------------------------
N_itration = 6
N_drugs = ncol(sen)
Results = data.frame()
#i=325
for (i in 1432:1433){
  print(paste0("The drug number is: ", as.character(i)))
  
  X = omic[!is.na(sen[,i]),]
  y = sen[!is.na(sen[,i]),i]
  
  # Cross validation loop
  
  clusterExport(cl, c("X","y","i","my_genes","MyGraph"))
  clusterEvalQ(cl, c(source("high_corr_FS.R"), source("mRMR.R"),
                     library(igraph),source("Infogenes_FS.R"),
                     library(caTools),library(randomForest),source("RF_Func_caret.R"),
                     library(glmnet),library(caret),source ("ENet_Func.R"),
                     library(keras),library(tensorflow),source("MLP_Func.R")))
  
  RepLoop = function(j){
    
    sample = sample.split(y, SplitRatio = .9)
    
    Xtr_val = subset(X, sample == TRUE)
    Xtest  = subset(X, sample == FALSE)
    ytr_val = subset(y, sample == TRUE)
    ytest  = subset(y, sample == FALSE)
    
    #sample = sample.split(ytr_val, SplitRatio = .8)
    
    #Xtrain = subset(Xtr_val, sample == TRUE)
    #Xval  = subset(Xtr_val, sample == FALSE)
    #ytrain = subset(ytr_val, sample == TRUE)
    #yval  = subset(ytr_val, sample == FALSE)
    
    Xtrain = Xtr_val
    ytrain = ytr_val
    
    # Normalization-------------------------------------------------------------
    # Xtrain normalization
    Mean_X = apply(Xtrain,2,mean)
    STD_X = apply(Xtrain,2,sd)
    Xtrain = (Xtrain-Mean_X)/STD_X
    
    # Xtest normalization
    Xtest = (Xtest-Mean_X)/STD_X
    
    # Ytrain normalization
    Mean_y = mean(ytrain)
    STD_y = sd(ytrain)
    ytrain_norm = (ytrain-Mean_y)/STD_y
    
    
    # Feature selection---------------------------------------------------------
    FS_method = "Infogenes"
      
    if(FS_method == "None"){
      Xtrain = Xtrain
      
    }else if(FS_method == "high_corr"){
      Xtrain = high_corr(Xtrain,ytrain,N_feat=100)
      
    }else if(FS_method == "mRMR"){
      Xtrain = mRMR(Xtrain, ytrain, N_feat = 1000, alpha=1, do.plot = FALSE)
      
    }else if(FS_method == "Infogenes"){
      Xtrain = Infogenes(Xtrain,ytrain,MyGraph,my_genes)
      
    }else if(FS_method == "DoRothEA"){
      
    }else if(FS_method == "Progeny"){
      
    }
    
    Xtest = Xtest[,colnames(Xtrain)]
    
    # Models
    
    y_pred_RF = RandomForest(ytrain = ytrain ,Xtrain = Xtrain,
                             Xtest = Xtest)
    
    y_pred_ENet = ElasticNet(ytrain = ytrain ,Xtrain = Xtrain,
                             Xtest = Xtest)
    
    y_pred_MLP = MLP(ytrain = ytrain ,Xtrain = Xtrain,
                     Xtest = Xtest)
    
    
    # y_pred re-normalization
    y_pred_RF = (y_pred_RF*STD_y)+Mean_y
    y_pred_ENet = (y_pred_ENet*STD_y)+Mean_y
    y_pred_MLP = (y_pred_MLP*STD_y)+Mean_y
    
    # Evaluation
    mse_RF = mean((ytest-y_pred_RF)^2)
    corr_RF = cor(ytest,y_pred_RF)
    
    mse_ENet = mean((ytest-y_pred_ENet)^2)
    corr_ENet = cor(ytest,y_pred_ENet)
    
    mse_MLP = mean((ytest-y_pred_MLP)^2)
    corr_MLP = cor(ytest,y_pred_MLP)
    
    
    
    result = data.frame(mse_RF = mse_RF,corr_RF = corr_RF,
                        mse_ENet = mse_ENet, corr_ENet = corr_ENet,
                        mse_MLP = mse_MLP, corr_MLP = corr_MLP)
    
    return(result)
    }

  
  result = parLapply(cl, sapply(1:N_itration, list), RepLoop) 
  
  Result = data.frame()
  for (k in 1:N_itration){
    Result = rbind(Result, result[[k]])
  }
  
  Results = rbind(Results, data.frame(mean_mse_RF = mean(Result$mse_RF),
                                      sd_mse_RF = sd(Result$mse_RF),
                                      mean_corr_RF = mean(Result$corr_RF),
                                      sd_corr_RF = sd(Result$corr_RF),
                                      mean_mse_ENet = mean(Result$mse_ENet),
                                      sd_mse_ENet = sd(Result$mse_ENet),
                                      mean_corr_ENet = mean(Result$corr_ENet),
                                      sd_corr_ENet = sd(Result$corr_ENet),
                                      mean_mse_MLP = mean(Result$mse_MLP),
                                      sd_mse_MLP = sd(Result$mse_MLP),
                                      mean_corr_MLP = mean(Result$corr_MLP),
                                      sd_corr_MLP = sd(Result$corr_MLP)))
  
}
stopCluster(cl)
Results = t(Results)
#colnames(Results) = paste0("drug_",c(1:2))

#saveRDS(Results,"Processed_Data/Result_All_Drugs.rds")
#cor_sort = data.frame(sort(Results[,3],decreasing = TRUE))

