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

Drug_Targets = readRDS("Processed_data/S1/drug_targets.rds")
rownames(Drug_Targets) = Drug_Targets$name  

Models = c("RandomForest","ElasticNet", "Lasso","Ridge","MLP")

clusterExport(cl, c("GE_PRISM","GE_TCGA","sen_PRISM","res_TCGA","Drug_Targets"))
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
  

  Drug_Targets_i = Drug_Targets[i,2]
  Drug_Targets_i = unlist(strsplit(Drug_Targets_i,", "))
  
  #if(length(Drug_Targets_i)!=0){
    if(length(Drug_Targets_i)>1){
      
      I_G = intersect(Drug_Targets_i,colnames(Xtest))
      Xtrain_reduced = as.matrix(Xtrain[,I_G])
      Xtest_reduced = as.matrix(Xtest[,I_G])
      #I_G = intersect(Drug_Targets_i,colnames(Xtest))
      #Xtrain_reduced = as.matrix(Xtrain[,I_G])
      #Xtest_reduced = as.matrix(Xtest[,I_G])
      
      #Xtrain_reduced = cbind(Xtrain_reduced,Xtrain_reduced)
      #Xtest_reduced = cbind(Xtest_reduced,Xtest_reduced)
    
    # }else if(length(Drug_Targets_i)>1){
    #   
    #   I_G = intersect(Drug_Targets_i,colnames(Xtest))
    #   Xtrain_reduced = as.matrix(Xtrain[,I_G])
    #   Xtest_reduced = as.matrix(Xtest[,I_G])
    # }
    
    # Ytrain normalization
    ytrain = scale(ytrain)
    ytrain = ytrain[,1]
    # Models
    y_pred = Ridge(ytrain = ytrain, Xtrain = Xtrain_reduced, Xtest = Xtest_reduced)
    
    # Evaluation
    pred = prediction(y_pred, ytest==1)
    AUC = performance(pred, measure = "auc")
    AUC = as.numeric(AUC@y.values)
    
    corr = cor(ytest,y_pred)
    ttest = t.test(y_pred[ytest==1], y_pred[ytest==2], alternative="greater")$p.value
    Ranksum = wilcox.test(y_pred[ytest==1], y_pred[ytest==2], alternative ="greater")$p.value
    
    result = data.frame(AUC = AUC, corr = corr, ttest=ttest, Ranksum = Ranksum)
  }else{
    result = data.frame(AUC = 0, corr = 0, ttest=1, Ranksum = 1)
  }
  
  return(result)
}

drugs = colnames(sen_PRISM)
result = parLapply(cl, sapply(drugs, list), DrugLoop) 

Result = data.frame()
for (k in drugs){
  Result = rbind(Result, result[[k]])
}

stopCluster(cl)

saveRDS(Result,"Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF9-DrugTargets_Ridge.rds")
print(sum(Result$Ranksum<0.05))
print(which(Result$Ranksum<0.05))
print(which(Result$ttest<0.05))

#boxplot(y_pred[ytest==1],y_pred[ytest==2])


