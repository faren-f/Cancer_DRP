rm(list=ls())
library(parallel)
library(caTools)
library(igraph)
library(dorothea)
source("F4-DoRothEA.R")
source("F9-decoupleR.R")

no_cores = detectCores()
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

#Read data--------------------------------------------------
#PRISM
#sen = readRDS("Processed_Data/S1/sensitivity_matrix_AUC.rds")
sen = readRDS("Processed_Data/Other/sen_reduced.rds")
GE = readRDS("Processed_Data/S1/expresion_matrix.rds")
drug_targets = readRDS("Processed_data/S1/drug_targets.rds")
#res_drugs = readRDS("Processed_data/Other/Result_All_Drugs.rds")
#order_drugs = data.frame(order = order(res_drugs[,3],decreasing = TRUE))
#drugs = data.frame(colnames(sen))

#CCLE
#sen = readRDS("Processed_Data/S15/sensitivity_matrix_Activity_Area.rds")
#GE = readRDS("Processed_Data/S15/expresion_matrix.rds")

Interaction_Network = "OP_decoupleR"
MyGraph = NA
my_genes = NA
if (Interaction_Network == "STRING"){
  
  GE = readRDS("Processed_data/S9/expresion_matrix_PRISM_STRING.rds")
  STRING_edgelist = readRDS("Processed_data/S9/ppi_STRING_PRISM.rds")
  ppi = rbind(STRING_edgelist$gene_symbol1, STRING_edgelist$gene_symbol2)
  MyGraph = simplify(graph(ppi, directed = FALSE))
  my_genes = V(MyGraph)$name
  omics = GE
  
} else if(Interaction_Network == "Omnipath"){
  
  GE = readRDS("Processed_data/S10/expresion_matrix_PRISM_Omnipath.rds")
  Omnipath_edgelist = readRDS("Processed_data/S10/ppi_Omnipath_PRISM.rds")
  ppi = rbind(Omnipath_edgelist$gene_symbol1, Omnipath_edgelist$gene_symbol2)
  MyGraph = simplify(graph(ppi, directed = FALSE))
  my_genes = V(MyGraph)$name
  omics = GE
  
} else if(Interaction_Network == "OP_decoupleR"){
  TF_dR = decoupleR(X = GE, method = "gsva")
  omics = TF_dR
  
} else if(Interaction_Network == "OP_DoRothEA"){
  TF_Do = DoRothEA(X = GE)
  omics = t(TF_Do)
}

  
#loop across drugs--------------------------------------
N_itration = 60
N_drugs = ncol(sen)
Results = c()

for (i in 1:N_drugs){
  print(paste0("The drug number is: ", as.character(i)))
  
  #drug_targets_i = strsplit(drug_targets[i,2],", ")
  #DT = intersect(colnames(GE),drug_targets_i[[1]])
  
  #if (length(DT)<1){
    #Results = c(Results,0)
    #next
  #}
  
  #GE_DT = GE[, c(DT, DT[1])]
  #omics = GE_DT
  
  X = omics[!is.na(sen[,i]),]
  y = sen[!is.na(sen[,i]),i]
  
  # Cross validation loop
  cl = makeCluster(no_cores-2)
  clusterExport(cl, c("X","y","i","my_genes","MyGraph","Interaction_Network"))
  clusterEvalQ(cl, c(source("F1-high_corr.R"), source("F2-mRMR.R"),
                     library(igraph),source("F3-Infogenes.R"),
                     library(dorothea),source("F4-DoRothEA.R"),source("F9-decoupleR.R"),
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
    FS_method = "OP_decoupleR"
    
    if (Interaction_Network == "STRING" | Interaction_Network=="Omnipath"){
      Xtrain = Infogenes(Xtrain,ytrain,MyGraph,my_genes)
      
    }else if(Interaction_Network == "OP_decoupleR"){
      Xtrain = Xtrain
      
    }else if(Interaction_Network == "" | Interaction_Network == "OP_DoRothEA"){
      
      if(FS_method == "high_corr"){
        Xtrain = high_corr(Xtrain,ytrain,N_feat=400)
        
      }else if(FS_method == "mRMR"){
        Xtrain = mRMR(Xtrain, ytrain, N_feat = 100, alpha=1, do.plot = FALSE)
        
      }else if(FS_method == "Progeny"){
        
      }
    }
    
    Xtest = Xtest[,colnames(Xtrain)]
    
    # Models
    
    y_pred_RF = RandomForest(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    
    #y_pred_ENet = Ridge(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    
    #y_pred_MLP = MLP(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    
    
    # y_pred re-normalization
    y_pred_RF = (y_pred_RF*STD_y)+Mean_y
    #y_pred_ENet = (y_pred_ENet*STD_y)+Mean_y
    #y_pred_MLP = (y_pred_MLP*STD_y)+Mean_y
    
    # Evaluation
    mse_RF = mean((ytest-y_pred_RF)^2)
    corr_RF = cor(ytest,y_pred_RF)
    
    #mse_ENet = mean((ytest-y_pred_ENet)^2)
    #corr_ENet = cor(ytest,y_pred_ENet)
    
    #mse_MLP = mean((ytest-y_pred_MLP)^2)
    #corr_MLP = cor(ytest,y_pred_MLP)
    
    
    result = data.frame(corr_RF = corr_RF, mse_RF = mse_RF)
                        #mse_ENet = mse_ENet, corr_ENet = corr_ENet)
                        #mse_MLP = mse_MLP, corr_MLP = corr_MLP)
    
    return(result)
  }
  
  
  result = parLapply(cl, sapply(1:N_itration, list), RepLoop) 
  
  Result = data.frame()
  for (k in 1:N_itration){
    Result = rbind(Result, result[[k]])
  }
  
  
  r = c(mean(Result$corr_RF))
                 #mean_corr_ENet = mean(Result$corr_ENet))
  print(r)
  Results = c(Results, r)
  
  
  #Results = rbind(Results, data.frame(mean_corr_RF = mean(Result$corr_RF),
                                      #mean_mse_RF = mean(Result$mse_RF),
                                      #sd_mse_RF = sd(Result$mse_RF),
                                      #sd_corr_RF = sd(Result$corr_RF),
                                      #mean_mse_ENet = mean(Result$mse_ENet),
                                      #sd_mse_ENet = sd(Result$mse_ENet),
                                      # mean_corr_ENet = mean(Result$corr_ENet)))
                                      #sd_corr_ENet = sd(Result$corr_ENet)))
                                      #mean_mse_MLP = mean(Result$mse_MLP),
                                      #sd_mse_MLP = sd(Result$mse_MLP),
                                      #mean_corr_MLP = mean(Result$corr_MLP),
                                      #sd_corr_MLP = sd(Result$corr_MLP)))
  

  stopCluster(cl)
}
#colnames(Results) = paste0("drug_",c(1:2))

#saveRDS(Results,"Processed_Data/Result_All_Drugs.rds")
#cor_sort = data.frame(sort(Results[,3],decreasing = TRUE))
Results = data.frame(Results)

####drug targets------------------------
#hist(Results[,1])
# Results =readRDS("Processed_data/S0/Results.rds")
# a = data.frame(order(Results[,1],decreasing = TRUE))
# b = data.frame(sort(Results[,1],decreasing = TRUE))
#"poziotinib"  "AZD8931"     "dacomitinib" "erdafitinib" "nutlin-3"    "BMS-690514" 
#"idasanutlin" "bosutinib"   "AMG-232"

#####------------------------------------------------------
rownames(Results) = 1:ncol(sen)
Results$drug_name = colnames(sen)
saveRDS(Results,"Processed_Data/S0/Result_decoupleR.rds")
e = readRDS("Processed_data/S0/Result_decoupleR.rds")
