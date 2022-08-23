rm(list=ls())
library(parallel)
library(igraph)
source("F3-Infogenes_DT.R")

no_cores = detectCores()
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

#Read data--------------------------------------------------
#Results_DT = readRDS("Processed_Data/S0/Results_drugTargets.rds")
sen = readRDS("Processed_Data/S0/sen_reduced.rds")

#sen = readRDS("Processed_Data/S1/sensitivity_matrix_AUC.rds")
#drugs = Results_DT[1:10,4]
#sen_re = sen[,intersect(colnames(sen),drugs)]

GE = readRDS("Processed_data/S9/expresion_matrix_PRISM_STRING.rds")
drug_targets = readRDS("Processed_data/S1/drug_targets.rds")
#drug_targets_re = drug_targets[drug_targets$name %in% drugs,]
drug_targets_re = drug_targets[drug_targets$name %in% colnames(sen),]
rownames(drug_targets_re) = drug_targets_re[,1]
drug_targets_re = drug_targets_re[colnames(sen),]

STRING_edgelist = readRDS("Processed_data/S9/ppi_STRING_PRISM.rds")
ppi = rbind(STRING_edgelist$gene_symbol1, STRING_edgelist$gene_symbol2)
MyGraph = simplify(graph(ppi, directed = FALSE))
my_genes = V(MyGraph)$name


#loop across drugs--------------------------------------
N_itration = 60
N_drugs = nrow(drug_targets_re)
Results = matrix(0,350,2)

for (i in 1:350){
  print(paste0("The drug number is: ", as.character(i)))
  
  drug_targets_i = strsplit(drug_targets_re[i,2],", ")
  DT = intersect(colnames(GE),drug_targets_i[[1]])

  if (length(DT)<1){
    Results[i,]= 0
    next
  }
  
  Score_diff = Infogenes(MyGraph,my_genes,DT = DT)
  
  X = GE[!is.na(sen[,i]),]
  y = sen[!is.na(sen[,i]),i]
  
  # Cross validation loop
  cl = makeCluster(no_cores-2)
  clusterExport(cl, c("X","y","i","Score_diff"))
  clusterEvalQ(cl, c(library(caTools),library(gelnet)))
  
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

    
    # Models
    d = rep(1, ncol(Xtrain))
    d_IG = 1-Score_diff
    n_feat = length(d)
    lambda1 = .5
    lambda2 = 0
    
    GelNet = gelnet(Xtrain, ytrain_norm, l1 = lambda1, l2 = lambda2, d = d,
                    P = diag(d), m = rep(0,n_feat), max.iter = 10, eps = 1e-05)
    
    GelNet_IG = gelnet(Xtrain, ytrain_norm, l1 = lambda1, l2 = lambda2, d = d_IG,
                              P = diag(d_IG), m = rep(0,n_feat), max.iter = 10, eps = 1e-05)
    
    Beta = GelNet[["w"]]
    Beta0 = GelNet[["b"]]
    y_pred = (Xtest %*%  Beta) + Beta0
    
    
    Beta_IG = GelNet_IG[["w"]]
    Beta0_IG = GelNet_IG[["b"]]
    y_pred_IG = (Xtest %*%  Beta_IG) + Beta0_IG
    
    
    # y_pred re-normalization
    y_pred = (y_pred*STD_y)+Mean_y
    y_pred_IG = (y_pred_IG*STD_y)+Mean_y
    
    
    # Evaluation
    mse = mean((ytest-y_pred)^2)
    corr = cor(ytest,y_pred)
    
    mse_IG = mean((ytest-y_pred_IG)^2)
    corr_IG = cor(ytest,y_pred_IG)
    result = data.frame(corr = corr, corr_IG = corr_IG, mse = mse, mse_IG = mse_IG)
    
    return(result)
  }
  
  result = parLapply(cl, sapply(1:N_itration, list), RepLoop) 
  
  Result = data.frame()
  for (k in 1:N_itration){
    Result = rbind(Result, result[[k]])
  }
  
  cor = c(mean(Result$corr))
  cor_IG = c(mean(Result$corr_IG))
  
  print(c(cor = cor,cor_IG = cor_IG))
  Results[i,] = c(cor, cor_IG)
  
  stopCluster(cl)
}

#saveRDS(Results,"Processed_Data/S0/Result_Novel_method.rds")

