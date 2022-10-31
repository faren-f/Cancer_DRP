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
if(sum(q3_genes==0)>0){
  GE_TCGA = GE_TCGA[,-which(q3_genes==0)]
  GE_PRISM = GE_PRISM[,-which(q3_genes==0)]
}

#Models = c("MLP")
Models = c("LinearcRegresion", "RandomForest","ElasticNet", "Lasso","Ridge","MLP")


clusterExport(cl, c("GE_PRISM","GE_TCGA","sen_PRISM","res_TCGA","Models"))
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
                   source("F19-Drug_Pathway_gene_set.R"),
                   source("F22-Drug_Pathway_Level_genes_eachTarget.R")))


DrugLoop = function(i){
  
  pathway_gene_set = Drug_Pathway_gene_set_eachTarget(drug = i, level=1)[["all"]]
  result = c()
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
      
      # Ytrain normalization
      ytrain = scale(ytrain)
      ytrain = ytrain[,1]
      
      # Models
      
      for(M in Models){
        model = get(M)
        y_pred = model(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
        corr = cor(ytest,y_pred)
        ttest = t.test(y_pred[ytest==1], y_pred[ytest==2], alternative="greater")$p.value
        Ranksum = wilcox.test(y_pred[ytest==1], y_pred[ytest==2], alternative ="greater")$p.value
        result = rbind(result, c(corr, ttest, Ranksum, N_genes))
      }
      
    }else{
      corr = 0
      ttest = 1
      Ranksum = 1
      N_genes = 0
      
      for(rep in 1:length(Models)){
        result = rbind(result, c(corr, ttest, Ranksum, N_genes))
      }
    }
  }else{
    corr = 0
    ttest = 1
    Ranksum = 1
    N_genes = 0
    for(rep in 1:length(Models)){
      result = rbind(result, c(corr, ttest, Ranksum, N_genes))
    }
  }
  
  return(result)
}


drugs = colnames(sen_PRISM)
result = parLapply(cl, sapply(drugs, list), DrugLoop) 


Result = data.frame()
for (k in drugs){
  Result = rbind(Result, result[[k]])
}

N_Models = length(Models)
Result_each_Model = list()
for(m in 1:N_Models){
  R = Result[m,]
  for(d in seq(N_Models,(N_Models*length(drugs))-1,N_Models)){
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

saveRDS(Result_each_Model, "Final_Result/Train@PRISM_Test@TCGA_Models/RF26-allModels_Pathway.rds")
# 5  8 16 17 24
#5 17
#3  5 24
