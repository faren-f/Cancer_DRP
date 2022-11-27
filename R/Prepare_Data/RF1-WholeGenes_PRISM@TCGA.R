rm(list=ls())

library(parallel)
no_cores = detectCores()
cl = makeCluster(no_cores-2)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
sen_PRISM = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")
res_TCGA = readRDS("Processed_data/Other/Res_TCGA_24_Drugs.rds")

GE_PRISM = readRDS("Processed_data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")
GE_TCGA = readRDS("Processed_data/S23/expresion_matrix_TCGA.rds")

# Remove genes whose Q3 is zero
q3_genes = apply(GE_TCGA,2,quantile,prob=0.75)
GE_TCGA = GE_TCGA[,-which(q3_genes==0)]
GE_PRISM = GE_PRISM[,-which(q3_genes==0)]

Models = c("RandomForest","ElasticNet", "Lasso","Ridge","MLP")

clusterExport(cl, c("GE_PRISM","GE_TCGA","sen_PRISM","res_TCGA"))
clusterEvalQ(cl, c(library(ROCR), source("F18-Combat_Normalization.R"),
                   source("F10-Ridge.R"),
                   source("F7-RandomForest.R"),
                   source("F6-ENet.R"),
                   source("F8-MLP.R"),
                   source("F13-Lasso.R")))

DrugLoop = function(i){
  
  print(paste0("The drug number is: ", as.character(i)))

  Xtrain = GE_PRISM[!is.na(sen_PRISM[,i]),]
  ytrain = sen_PRISM[!is.na(sen_PRISM[,i]),i]
  
  Xtest = GE_TCGA[!is.na(res_TCGA[,i]),]
  ytest = res_TCGA[!is.na(res_TCGA[,i]),i]
  
  X_Normalization = Combat_Scale(Xtrain,Xtest)
  
  Xtrain = X_Normalization[[1]]
  Xtest = X_Normalization[[2]]
  
  source("F15-Feature_Selection_PRISM@TCGA.R")
  selected_features = c("Whole_genes")
  Omics_List = Feature_Selection_PRISM_TCGA(selected_features, Xtrain=Xtrain, Xtest=Xtest)
  Xtrain = Omics_List[[1]]
  index = Omics_List[[2]]
  Xtest = Omics_List[[3]]
  
  # Ytrain normalization
  ytrain = scale(ytrain)
  ytrain = ytrain[,1]
  # Models
  y_pred = MLP(ytrain = ytrain ,Xtrain = Xtrain, Xtest = Xtest)

  # Evaluation
  pred = prediction(y_pred, ytest==1)
  AUC = performance(pred, measure = "auc")
  AUC = as.numeric(AUC@y.values)
  
  corr = cor(ytest,y_pred)
  ttest = t.test(y_pred[ytest==1], y_pred[ytest==2], alternative="greater")$p.value
  Ranksum = wilcox.test(y_pred[ytest==1], y_pred[ytest==2], alternative ="greater")$p.value
  
  result = data.frame(AUC = AUC, corr = corr, ttest=ttest, Ranksum = Ranksum)
  
  return(result)
}

N_drug = ncol(sen_PRISM)
result = parLapply(cl, sapply(1:N_drug, list), DrugLoop) 

Result = data.frame()
for (k in 1:N_drug){
  Result = rbind(Result, result[[k]])
}

stopCluster(cl)

saveRDS(Result,"Final_Result/TrainPRISM&TestTCGA_FS/MLP/RF1_WholeGenes_MLP.rds")
print(sum(Result$Ranksum<0.05))
print(which(Result$Ranksum<0.05))
print(which(Result$ttest<0.05))

