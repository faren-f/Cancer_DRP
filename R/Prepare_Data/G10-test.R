rm(list=ls())

library(parallel)
no_cores = detectCores()
cl = makeCluster(no_cores-2)

source("F19-Drug_Pathway_gene_set.R")
source("F14-Feature_Selection.R")

GE = readRDS("Processed_Data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")
TCGA_good_drugs = c("bicalutamide", "docetaxel", "etoposide", "paclitaxel", "leucovorin", 
                    "dacarbazine", "methotrexate", "ifosfamide", "gemcitabine", 
                    "vincristine", "cisplatin","vinblastine")
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
#sen = readRDS("All_Results/sen_PRISM_good_drugs.rds")
sen = readRDS("Processed_data/S1/sensitivity_matrix_AUC.rds")
#TCGA_PRISM_drugs_all = readRDS("Processed_data/S21/Drugs_TCGA@PRISM.rds")
TCGA_PRISM_drugs_sig_samples = readRDS("Processed_data/Other/PRISM_TCGA_drugs.rds")

which(colnames(sen) %in% TCGA_PRISM_drugs_sig_samples)
I =intersect(colnames(sen),TCGA_PRISM_drugs_sig_samples)
sen = sen[,I]

#Feature selection
# selected_features = c("Landmark_genes")
# Omics_List = Feature_Selection(selected_features,GE)
# X = Omics_List[[1]]
# index = Omics_List[[2]]

drugs = data.frame(colnames(sen))
N_drug = ncol(sen)
Results = c()
N_genes = c()
C = c()

#4,5,14,16,17,31,39,48,75,78
for (i in 21){
  print(paste0("The drug number is: ", as.character(i)))
  drug = drugs[i,1]
  
  #Drug Pathway feature selection
  pathway_gene_set = Drug_Pathway_gene_set(drug)
  
  if(isEmpty(pathway_gene_set[[1]])){
    C = c(C,i)
    print(i)
    next
    }
    
  I = intersect(colnames(GE),pathway_gene_set[,1])
  X = GE[,I]
  N_genes = c(N_genes, length(I))
  index = rep(1,ncol(X))

  X = X[!is.na(sen[,i]),]
  y = sen[!is.na(sen[,i]),i]

  #X = scale(X)

  # Ytrain normalization
  #Mean_y = mean(y)
  #STD_y = sd(y)
  #y = (y-Mean_y)/STD_y


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
    #corr_MLP = cor(ytest,y_pred_MLP)
      
    result = data.frame(corr_Ridge = corr_Ridge)
      #corr_SGL = corr_SGL,
      #corr_RF = corr_RF,
      #corr_ENet = corr_ENet,
      #corr_Lasso = corr_Lasso,
      #corr_MLP = corr_MLP)
      
      return(result)
    }
    
    N_itration = 100
    result = parLapply(cl, sapply(1:N_itration, list), RepLoop) 
  
    Result = data.frame()
    for (k in 1:N_itration){
      Result = rbind(Result, result[[k]])
    }
  
    Result_mean = apply(Result, 2, mean)
    Result_sd = apply(Result, 2, sd)
    print(Result_mean)
    print(Result_sd)
    print(dim(X))
    
    Results = rbind(Results, c(Result_mean, Result_sd))
  
}
stopCluster(cl)
#c = colnames(sen)
#rownames(Results) = c[order_drugs[3:20,]]
#saveRDS(Results,"All_Results/.rds")
#all = readRDS("All_Results/SGL_RF@L1000_TF@12run.rds")
#a = data.frame(colnames(sen))


#s = colnames(sen)[which(Results[,1]>0.2)]
#intersect(s, TCGA_good_drugs)

