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

N_Level = 10 
drugs = colnames(sen_PRISM)

Models = c("RandomForest","ElasticNet", "Lasso","Ridge","MLP")

clusterExport(cl, c("GE_PRISM","GE_TCGA","sen_PRISM","res_TCGA","drugs","N_Level"))
clusterEvalQ(cl, c(library(ROCR), source("F18-Combat_Normalization.R"),
                   source("F10-Ridge.R"),
                   source("F7-RandomForest.R"),
                   source("F6-ENet.R"),
                   source("F8-MLP.R"),
                   source("F13-Lasso.R")))

DrugLoop = function(i){
  
  #Drug Pathway feature selection
  Result_Level = list()
  for (j in 1:N_Level){
    
    source("F22-Drug_Pathway_Level_genes_eachTarget.R")
    pathway_gene_set = Drug_Pathway_gene_set_eachTarget(drug = i, level=j)
    
    if(!isEmpty(pathway_gene_set)){
       
      Result = c()
      for(k in names(pathway_gene_set)){
        
        if(!is.null(pathway_gene_set[[k]])){
          I = intersect(colnames(GE_PRISM),pathway_gene_set[[k]])
          X_PRISM = GE_PRISM[,I]
          X_TCGA = GE_TCGA[,I]
          if(!is.null(ncol(X_PRISM))){
            
            Xtrain = X_PRISM[!is.na(sen_PRISM[,i]),]
            ytrain = sen_PRISM[!is.na(sen_PRISM[,i]),i]
            
            Xtest = X_TCGA[!is.na(res_TCGA[,i]),]
            ytest = res_TCGA[!is.na(res_TCGA[,i]),i]
            
            X_Normalization = Combat_Scale(Xtrain,Xtest)
            
            Xtrain = X_Normalization[[1]]
            Xtest = X_Normalization[[2]]
            N_genes = ncol(Xtrain)
            
            # Ytrain normalization
            ytrain = scale(ytrain)
            ytrain = ytrain[,1]
    
            # Models
            y_pred = RandomForest(ytrain = ytrain ,Xtrain = Xtrain, Xtest = Xtest)
            
            # Evaluation
            pred = prediction(y_pred, ytest==1)
            AUC = performance(pred, measure = "auc")
            AUC = as.numeric(AUC@y.values)
            
            corr = cor(ytest,y_pred)
            ttest = t.test(y_pred[ytest==1], y_pred[ytest==2], alternative="greater")$p.value
            Ranksum = wilcox.test(y_pred[ytest==1], y_pred[ytest==2], alternative ="greater")$p.value
            Result = rbind(Result, cbind(AUC, corr, ttest, Ranksum, N_genes))
            
          }else{
            AUC = 0
            N_genes = 0
            corr = 0
            ttest = 1
            Ranksum = 1
            Result = rbind(Result, cbind(AUC, corr,ttest,Ranksum,N_genes))
          }
        }else{
          AUC = 0
          N_genes = 0
          corr = 0
          ttest = 1
          Ranksum = 1
          Result = rbind(Result, cbind(AUC, corr,ttest,Ranksum,N_genes))
        }
      }
      rownames(Result) = names(pathway_gene_set)
    }else{
      Result = c()
    }
    Result_Level[j] = list(Result)
  }
  Result_Drug = Result_Level
  #result = cbind(Result, Targets = names(pathway_gene_set))
  return(Result_Drug)
}

Result_Drug = parLapply(cl, sapply(drugs, list), DrugLoop) 

stopCluster(cl)

saveRDS(Result_Drug,"Final_Result/TrainPRISM&TestTCGA_FS/RF/RF6-PW_EachTarget_RF.rds")


