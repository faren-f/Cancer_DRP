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
clusterEvalQ(cl, c(source("F18-Combat_Normalization.R"),
                   source("F10-Ridge.R"),
                   source("F7-RandomForest.R"),
                   source("F6-ENet.R"),
                   source("F8-MLP.R"),
                   source("F13-Lasso.R"),
                   source("F22-Drug_Pathway_Level_genes_eachTarget.R")))
             
DrugLoop = function(i){
  
  pathway_gene_set = Drug_Pathway_gene_set_eachTarget(drug = i, level=1)[["all"]]
  
  if(!isEmpty(pathway_gene_set)){

    I = intersect(colnames(GE_PRISM),pathway_gene_set)
    X = GE_PRISM[,I]
    X_TCGA = GE_TCGA[,I]
    N_genes = length(I)

    Xtrain = X[!is.na(sen_PRISM[,i]),]
    ytrain = sen_PRISM[!is.na(sen_PRISM[,i]),i]
    
    Xtest = X_TCGA[!is.na(res_TCGA[,i]),]
    ytest = res_TCGA[!is.na(res_TCGA[,i]),i]
    
    length(ytest)
    if(length(ytest)>10){
      
      X_Normalization = Combat_Scale(Xtrain,Xtest)
      
      Xtrain = X_Normalization[[1]]
      Xtest = X_Normalization[[2]]
      
      y_pred = MLP(ytrain = ytrain ,Xtrain = Xtrain, Xtest = Xtest)
      corr = cor(ytest,y_pred)
      
      
      ttest = t.test(y_pred[ytest==1], y_pred[ytest==2], alternative="greater")$p.value
      Ranksum = wilcox.test(y_pred[ytest==1], y_pred[ytest==2], alternative ="greater")$p.value
  }else{
    corr = 0
    ttest = 1
    Ranksum = 1
    N_genes = 0
    }
  }else{
    corr = 0
    ttest = 1
    Ranksum = 1
    N_genes = 0
  }

  result = data.frame(corr = corr, ttest=ttest, 
                      Ranksum = Ranksum, N_genes = N_genes)
  
  return(result)
}

drugs = colnames(sen_PRISM)
result = parLapply(cl, sapply(drugs, list), DrugLoop) 

Result = c()
for (k in drugs){
  Result = rbind(Result, result[[k]])
}

stopCluster(cl)

saveRDS(Result,"Final_Result/TrainPRISM&TestTCGA_FS/MLP/RF6-PW_MLP.rds")
print(sum(Result$Ranksum<0.05))
print(which(Result$Ranksum<0.05))
print(which(Result$ttest<0.05))


