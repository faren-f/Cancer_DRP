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

clusterExport(cl, c("GE_PRISM","GE_TCGA","sen_PRISM","res_TCGA"))
clusterEvalQ(cl, c(source("F10-Ridge.R"), source("F18-Combat_Normalization.R"), source("F8-MLP.R")))

RepLoop = function(j){
  i = 16
  RandomSelectedGenes = sample(ncol(GE_PRISM),280)
  GE_PRISM_reduced = GE_PRISM[,RandomSelectedGenes]
  GE_TCGA_reduced = GE_TCGA[,RandomSelectedGenes]
  
  Xtrain = GE_PRISM_reduced[!is.na(sen_PRISM[,i]),]
  ytrain = sen_PRISM[!is.na(sen_PRISM[,i]),i]
  
  Xtest = GE_TCGA_reduced[!is.na(res_TCGA[,i]),]
  ytest = res_TCGA[!is.na(res_TCGA[,i]),i]
  
  X_Normalization = Combat_Scale(Xtrain,Xtest)
  Xtrain = X_Normalization[[1]]
  Xtest = X_Normalization[[2]]
  
  # Ytrain normalization
  ytrain = scale(ytrain)
  ytrain = ytrain[,1]
  # Models
  y_pred = Ridge(ytrain = ytrain, Xtrain = Xtrain, Xtest = Xtest)
  #y_pred = MLP(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
  
  # Evaluation
  corr = cor(ytest,y_pred)
  ttest = t.test(y_pred[ytest==1], y_pred[ytest==2], alternative="greater")$p.value
  Ranksum = wilcox.test(y_pred[ytest==1], y_pred[ytest==2], alternative ="greater")$p.value
  
  result = data.frame(corr = corr, ttest=ttest, Ranksum = Ranksum)
  
  return(result)
}

Rep = 1000
result = parLapply(cl, sapply(1:Rep, list), RepLoop) 

Result = data.frame()
for (k in 1:Rep){
  Result = rbind(Result, result[[k]])
}

stopCluster(cl)

#saveRDS(Result,"Final_Result/Train@PRISM_Test@TCGA_FS/RF0_Random.rds")
print(sum(Result$Ranksum<0.05))
print(which(Result$Ranksum<0.05))
print(which(Result$ttest<0.05))

