rm(list=ls())

library(parallel)
no_cores = detectCores()
cl = makeCluster(no_cores-2)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
sen_PRISM = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")
res_TCGA = readRDS("Processed_data/Other/Res_TCGA_24_Drugs.rds")

dR_PRISM = read.table("Processed_data/S33/gsea2_PRISM.csv",sep = ",",header = TRUE, row.names = 1)
dR_TCGA = read.table("Processed_data/S33/gsea2_TCGA.csv",sep = ",",header = TRUE, row.names = 1)

# Remove genes whose Q3 is zero
q3_genes = apply(dR_TCGA,2,quantile,prob=0.75)
if(sum(q3_genes==0)>0){
  dR_TCGA = dR_TCGA[,-which(q3_genes==0)]
  dR_PRISM = dR_PRISM[,-which(q3_genes==0)]
}


#Models = c("MLP")
Models = c("LinearcRegresion", "RandomForest","ElasticNet", "Lasso","Ridge","MLP")

clusterExport(cl, c("dR_PRISM","dR_TCGA","sen_PRISM","res_TCGA","Models"))
clusterEvalQ(cl, c(source("F7-RandomForest.R"),
                   source("F6-ENet.R"),
                   source("F8-MLP.R"),
                   source("F10-Ridge.R"),
                   source("F11-SGL.R"),
                   source("F13-Lasso.R"), 
                   source("F25-LinearRegression.R"),
                   source("F16-Zscore_Normalization.R"),
                   source("F17-Rank_Normalization.R"),
                   source("F18-Combat_Normalization.R"),
                   source("F19-Drug_Pathway_gene_set.R")))

DrugLoop = function(i){
  
  Xtrain = dR_PRISM[!is.na(sen_PRISM[,i]),]
  ytrain = sen_PRISM[!is.na(sen_PRISM[,i]),i]
  
  Xtest = dR_TCGA[!is.na(res_TCGA[,i]),]
  ytest = res_TCGA[!is.na(res_TCGA[,i]),i]
  
  #X_Normalization = Rank(Xtrain,Xtest)
  #X_Normalization = Rank(Xtrain,Xtest)
  X_Normalization = Combat_Scale(Xtrain,Xtest)
  
  Xtrain = X_Normalization[[1]]
  Xtest = X_Normalization[[2]]
  
  # Ytrain normalization
  ytrain = scale(ytrain)
  ytrain = ytrain[,1]
  
  # Models
  result = c()
  for(M in Models){
    model = get(M)
    y_pred = model(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    corr = cor(ytest,y_pred)
    ttest = t.test(y_pred[ytest==1], y_pred[ytest==2], alternative="greater")$p.value
    Ranksum = wilcox.test(y_pred[ytest==1], y_pred[ytest==2], alternative ="greater")$p.value
    result = rbind(result, c(corr, ttest, Ranksum))
  }
  
  return(result)
}

N_drug = ncol(sen_PRISM)
result = parLapply(cl, sapply(1:N_drug, list), DrugLoop) 


#R = unlist(result)

Result = data.frame()
for (k in 1:N_drug){
  Result = rbind(Result, result[[k]])
}

N_Models = length(Models)
Result_each_Model = list()
for(m in 1:N_Models){
  R = Result[m,]
  for(d in seq(N_Models,(N_Models*N_drug)-1,N_Models)){
    R = rbind(R, Result[m+d,])
  }
  Result_each_Model[m] = list(R)
}

stopCluster(cl)
Sum_Sig = c()
for(p in 1:N_Models){
  
  Sum_Sig = c(Sum_Sig, sum(Result_each_Model[[p]][3]<0.05))
  print(which(Result_each_Model[[p]][3]<0.05))
  #print(which(Result_each_Model[[p]][2]<0.05))
}

saveRDS(Result_each_Model, "Final_Result/Train@PRISM_Test@TCGA_Models/RF25-allModels_dR@gsea.rds")


