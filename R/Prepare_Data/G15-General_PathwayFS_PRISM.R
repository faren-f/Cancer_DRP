rm(list=ls())
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
source("F10-Ridge.R")
source("F7-RandomForest.R")
source("F8-MLP.R")
source("F13-Lasso.R")
source("F6-ENet.R")

#sen = readRDS("All_Results/sen_PRISM_good_drugs.rds")
GE = readRDS("Processed_Data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")
sen = readRDS("Processed_data/S1/sensitivity_matrix_AUC.rds")
#TCGA_PRISM_drugs_all = readRDS("Processed_data/S21/Drugs_TCGA@PRISM.rds")
TCGA_PRISM_drugs_sig_samples = readRDS("Processed_data/Other/PRISM_TCGA_drugs.rds")

which(colnames(sen) %in% TCGA_PRISM_drugs_sig_samples)
I =intersect(colnames(sen),TCGA_PRISM_drugs_sig_samples)
sen = sen[,I]
drugs = data.frame(colnames(sen))

N_drug = ncol(sen)
Results = c()

for (i in 16){
  print(paste0("The drug number is: ", as.character(i)))
  drug = drugs[i,1]
  
  for(j in 5:5){
    #Drug Pathway feature selection
    source("F22-Drug_Pathway_Level_genes_eachTarget.R")
    pathway_gene_set = Drug_Pathway_gene_set_eachTarget(drug = drug, level=j)
    
    if(!isEmpty(pathway_gene_set)){
      
      Result = c()
      for(k in names(pathway_gene_set)){
        
        if(!is.null(pathway_gene_set[[k]])){
          I = intersect(colnames(GE), pathway_gene_set[[k]])
          X = GE[,I]
          if(!is.null(ncol(X))){
            index = rep(1,ncol(X))
            X = X[!is.na(sen[,i]),]
            y = sen[!is.na(sen[,i]),i]
            
            #X = t(scale(t(X)))
            #X = (X-min(X))/(max(X)-min(X))
            
            # Ytrain normalization
            #y = (y-min(y))/(max(y)-min(y))
            hist(y)
            sample = sample.split(y, SplitRatio = .8)
            
            Xtrain = subset(X, sample == TRUE)
            Xtest  = subset(X, sample == FALSE)
            ytrain = subset(y, sample == TRUE)
            ytest  = subset(y, sample == FALSE)
            
            # a = hist(ytrain)
            # InvCount = 1/a[["counts"]]
            # 
            # W = rep(0,length(ytrain))
            # b = round(a$breaks,1)
            # for(l in 1:length(ytrain)){
            #   for(u in 2:length(b)){
            #     if((b[u-1] < ytrain[l]) & (ytrain[l] < b[u])){
            #       W[l] = InvCount[u-1]
            #     }
            #   }
            # }
            # 
            
            # Models
            #y_pred_SGL = My_SGL(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest,index = index)
            #y_pred_Ridge = RandomForest(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
            y_pred_ENet = ElasticNet(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
            y_pred_Lasso = Lasso(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
            #y_pred_Ridge = Ridge(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
            y_pred_MLP = MLP(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
            
            # Evaluation
            #corr_SGL = cor(ytest,y_pred_SGL)
            #corr_RF = cor(ytest,y_pred_RF)
            #corr_ENet = cor(ytest,y_pred_ENet)
            #corr_Lasso = cor(ytest,y_pred_Lasso)
            corr_Ridge = cor(ytest,y_pred_Ridge,method = "pearson")
            corr_Ridge
            #corr_MLP = cor(ytest,y_pred_MLP)
            hist(ytest)
            hist(y_pred_Ridge)
            plot(ytest,y_pred_Ridge,xlim = c(0,1.5), ylim = c(0,1.5))
            
          }else{
            Result = 0
          }
        }else{
          Result = 0
        }
      }
    }else{
    Result = 0
    }
  }
}
  
 
    
    