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
res = readRDS("Processed_data/S0/Result_All_Drugs.rds")
res = order(res[,3],decreasing = TRUE)
sen = readRDS("Processed_Data/S1/sensitivity_matrix.rds")
GE = readRDS("Processed_Data/S1/expresion_matrix.rds")
##TF
#TF_Do = DoRothEA(X = GE)
#TF_Do = t(TF_Do)
#omics = TF_Do

#decoupleR
TF_dR = decoupleR(X = GE, method = "gsva")
omics = TF_dR


##L1000
l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
#GE = GE[,colnames(GE)%in%l1000_genes]

#sample_tissue_types = readRDS("Processed_data/S19/sample_tissue_types.rds")
#GE = cbind(TF_Do,GE)
#omics = GE

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
N_itration = 60
N_drugs = ncol(sen)
Results = data.frame()
#i=325
drug = 1253
for (i in drug:drug){
  print(paste0("The drug number is: ", as.character(i)))
  
  X = omics[!is.na(sen[,i]),]
  y = sen[!is.na(sen[,i]),i]
  
  # Cross validation loop
  
  clusterExport(cl, c("X","y","i","my_genes","MyGraph","Interaction_Network"))
  clusterEvalQ(cl, c(source("F1-high_corr.R"), source("F2-mRMR.R"),
                     library(igraph),source("F3-Infogenes.R"),
                     library(dorothea),source("F4-DoRothEA.R"),source("F9-decoupleR.R"),
                     library(caTools),library(randomForest),source("F7-RandomForest.R")))
  
  RepLoop = function(j){
    
    sample = sample.split(y, SplitRatio = .9)
    
    Xtrain = subset(X, sample == TRUE)
    Xtest  = subset(X, sample == FALSE)
    ytrain = subset(y, sample == TRUE)
    ytest  = subset(y, sample == FALSE)
    
    # Normalization-------------------------------------------------------------
    # Xtrain normalization
    #Mean_X = apply(Xtrain,2,mean)
    #STD_X = apply(Xtrain,2,sd)
    #Xtrain = (Xtrain-Mean_X)/STD_X
    
    # Xtest normalization
    #Xtest = (Xtest-Mean_X)/STD_X
    
    # Ytrain normalization
    #Mean_y = mean(ytrain)
    #STD_y = sd(ytrain)
    #ytrain_norm = (ytrain-Mean_y)/STD_y
    
    # Feature selection---------------------------------------------------------
    FS_method = ""
    
    if (Interaction_Network == "STRING" | Interaction_Network=="Omnipath"){
      Xtrain = Infogenes(Xtrain,ytrain,MyGraph,my_genes)
      
    }else if(Interaction_Network == "OP_decoupleR"){
      Xtrain = Xtrain
      
    }else if(Interaction_Network == "None" | Interaction_Network == "OP_DoRothEA"){
      
      if(FS_method == "high_corr"){
        Xtrain = high_corr(Xtrain,ytrain,N_feat=200)
        
      }else if(FS_method == "mRMR"){
        Xtrain = mRMR(Xtrain, ytrain, N_feat = 200, alpha=1, do.plot = FALSE)
        
      }
    }
    
    Xtest = Xtest[,colnames(Xtrain)]
    
    # Models
    
    y_pred_RF = RandomForest(ytrain = ytrain ,Xtrain = Xtrain,
                             Xtest = Xtest)
  
    
    # y_pred re-normalization
    #y_pred_RF = (y_pred_RF*STD_y)+Mean_y
     
    # Evaluation
    mse_RF = mean((ytest-y_pred_RF)^2)
    corr_RF = cor(ytest,y_pred_RF)
    
    
    result = data.frame(mse_RF = mse_RF,corr_RF = corr_RF)
    
    return(result)
  }
  
  result = parLapply(cl, sapply(1:N_itration, list), RepLoop) 
  
  Result = data.frame()
  for (k in 1:N_itration){
    Result = rbind(Result, result[[k]])
  }
  
  Results = rbind(Results, data.frame(mean_corr_RF = mean(Result$corr_RF),
                                      mean_mse_RF = mean(Result$mse_RF),
                                      sd_mse_RF = sd(Result$mse_RF),
                                      sd_corr_RF = sd(Result$corr_RF)))
  
}
stopCluster(cl)
Results = t(Results)

