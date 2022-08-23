rm(list=ls())
library(parallel)
library(caTools)
library(igraph)
library(dorothea)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
res = readRDS("Processed_data/Other/Result_All_Drugs.rds")
r = order(res[,3],decreasing = TRUE)
r = data.frame(r)

source("F4-DoRothEA.R")
source("F9-decoupleR.R")

no_cores = detectCores()
cl = makeCluster(no_cores-2)


#Read data--------------------------------------------------
sen = readRDS("Processed_Data/S1/sensitivity_matrix.rds")
GE = readRDS("Processed_Data/S1/expresion_matrix.rds")
omics = GE
drugs = data.frame(colnames(sen))


Interaction_Network = ""
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
  TF_dR = decoupleR(X = t(GE), method = "gsva")
  omics = TF_dR
  
} else if(Interaction_Network == "OP_DoRothEA"){
  TF_Do = DoRothEA(X = GE)
  omics = t(TF_Do)
}


#loop across drugs--------------------------------------
N_itration = 6
N_drugs = ncol(sen)
Results = data.frame()

drug = 533
for (i in drug:drug){
  print(paste0("The drug number is: ", as.character(i)))
  
  X = omics[!is.na(sen[,i]),]
  y = sen[!is.na(sen[,i]),i]
  
  # Cross validation loop
  
  clusterExport(cl, c("X","y","i","my_genes","MyGraph","Interaction_Network"))
  clusterEvalQ(cl, c(source("F1-high_corr.R"), source("F2-mRMR.R"),
                     library(igraph),source("F3-Infogenes.R"),
                     library(dorothea),source("F4-DoRothEA.R"),source("F9-decoupleR.R"),
                     library(caTools),library(randomForest),source("F7-RandomForest.R"),
                     library(keras),library(tensorflow),source("F8-MLP.R")))
  
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
    FS_method = c("mRMR","high_corr")
    
    if (Interaction_Network == "STRING" | Interaction_Network=="Omnipath"){
      Xtrain = Infogenes(Xtrain,ytrain,MyGraph,my_genes)
      
    }else if(Interaction_Network == "OP_decoupleR"){
      Xtrain = Xtrain
      
    }else if(Interaction_Network == "" | Interaction_Network == "OP_DoRothEA"){
      
      for(f in FS_method){
        
        if(f == "high_corr"){
          Xtr = high_corr(Xtrain,ytrain,N_feat=200)
          Xte = Xtest[,colnames(Xtr)]
          y_pred_RF = MLP(ytrain = ytrain ,Xtrain = Xtr,Xtest = Xte)
          y_pred_RF = (y_pred_RF*STD_y)+Mean_y
          mse_RF_hc = mean((ytest-y_pred_RF)^2)
          corr_RF_hc = cor(ytest,y_pred_RF)
          
        }else if(f == "mRMR"){
          Xtr = mRMR(Xtrain, ytrain, N_feat = 200, alpha=1, do.plot = FALSE)
          Xte = Xtest[,colnames(Xtr)]
          y_pred_RF = MLP(ytrain = ytrain ,Xtrain = Xtr,Xtest = Xte)
          y_pred_RF = (y_pred_RF*STD_y)+Mean_y
          mse_RF_mRMR = mean((ytest-y_pred_RF)^2)
          corr_RF_mRMR = cor(ytest,y_pred_RF)
          
        }
      }
    }
    
    #Xtest = Xtest[,colnames(Xtrain)]
    
    # Models
    
    #y_pred_RF = RandomForest(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    
    
    # y_pred re-normalization
    #y_pred_RF = (y_pred_RF*STD_y)+Mean_y
    
    # Evaluation
    #mse_RF = mean((ytest-y_pred_RF)^2)
    #corr_RF = cor(ytest,y_pred_RF)
    
    
    result = data.frame(mse_RF_hc = mse_RF_hc,mse_RF_mRMR = mse_RF_mRMR,
                        corr_RF_hc = corr_RF_hc, corr_RF_mRMR = corr_RF_mRMR)
    
    return(result)
  }
  
  result = parLapply(cl, sapply(1:N_itration, list), RepLoop) 
  
  Result = data.frame()
  for (k in 1:N_itration){
    Result = rbind(Result, result[[k]])
  }
  
  Results = rbind(Results, data.frame(mean_corr_RF_hc = mean(Result$corr_RF_hc),
                                      mean_corr_RF_mRMR = mean(Result$corr_RF_mRMR),
                                      mean_mse_RF_hc = mean(Result$mse_RF_hc),
                                      mean_mse_RF_mRMR = mean(Result$mse_RF_mRMR)
                                      ))
  
}
stopCluster(cl)
Results = t(Results)

