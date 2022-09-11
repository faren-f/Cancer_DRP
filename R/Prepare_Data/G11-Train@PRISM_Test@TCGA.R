rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

source("F7-RandomForest.R")
source("F6-ENet.R")
source("F8-MLP.R")
source("F10-Ridge.R")
source("F11-SGL.R")
source("F13-Lasso.R")
source("F16-Zscore_Normalization.R")
source("F17-Rank_Normalization.R")
source("F18-Combat_Normalization.R")

low_sample_drugs = c(1,3,4,5,6,7,12,13,14,19,20,21,22,23,24,25,28,29,31,33,36,37
                     ,39,40,41,42,43,45,47,49,51,52,56,57,58)
sen_PRISM = readRDS("Processed_data/S23/sensitivity_matrix_PRISM_with@TCGA@drugs.rds")
res_TCGA = readRDS("Processed_data/S24/Drug_response_TCGA_binarized.rds")
#res_TCGA = readRDS("Processed_data/S23/Drug_response_matrix_TCGA.rds")

sen_PRISM = sen_PRISM[,-low_sample_drugs]
res_TCGA = res_TCGA[,-low_sample_drugs]

GE = readRDS("Processed_data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")
GE_TCGA = readRDS("Processed_data/S23/expresion_matrix_TCGA.rds")

# Remove genes whose Q3 is zero
q3_genes = apply(GE_TCGA,2,quantile,prob=0.75)
GE_TCGA = GE_TCGA[,-which(q3_genes==0)]
GE = GE[,-which(q3_genes==0)]

N_drug = ncol(sen_PRISM)
drugs = data.frame(colnames(sen_PRISM))
#saveRDS(drugs,"Processed_data/Other/24_drugs.rds")
Results = c()
N_genes = c()
C = c()
for (i in 1:N_drug){
  print(paste0("The drug number is: ", as.character(i)))
  drug = drugs[i,1]
  
  #Drug Pathway feature selection
  source("F19-Drug_Pathway_gene_set.R")
  pathway_gene_set = Drug_Pathway_gene_set(drug)

  if(isEmpty(pathway_gene_set[[1]])){
    print(i)
    C = c(C,i)
    next
    }
    
  I = intersect(colnames(GE),pathway_gene_set[,1])
  X = GE[,I]
  X_TCGA = GE_TCGA[,I]
  N_genes = c(N_genes, length(I))
  index = rep(1,ncol(X))
  ##################
  
  Xtrain = X[!is.na(sen_PRISM[,i]),]
  ytrain = sen_PRISM[!is.na(sen_PRISM[,i]),i]
  
  Xtest = X_TCGA[!is.na(res_TCGA[,i]),]
  ytest = res_TCGA[!is.na(res_TCGA[,i]),i]
  
  length(ytest)
  if(length(ytest)>10){
    
    #X_Normalization = Rank(Xtrain,Xtest)
    #X_Normalization = Rank(Xtrain,Xtest)
    X_Normalization = Combat_Scale(Xtrain,Xtest)
    
    Xtrain = X_Normalization[[1]]
    Xtest = X_Normalization[[2]]
    
    #source("F15-Feature_Selection_PRISM@TCGA.R")
    #selected_features = c("Whole_genes")
    #Omics_List = Feature_Selection(selected_features,GE = Xtrain ,GE_test = Xtest)
    # Xtrain = Omics_List[[1]]
    # index = Omics_List[[2]]
    # Xtest = Omics_List[[3]]
    
    # Ytrain normalization
    # Mean_ytrain = mean(ytrain)
    # STD_ytrain = sd(ytrain)
    # ytrain = (ytrain-Mean_ytrain)/STD_ytrain

    # Models
    #y_pred_Ridge = My_SGL(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest,index = index)
    #y_pred_Ridge = RandomForest(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    y_pred_Ridge = ElasticNet(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    #y_pred_Ridge = Lasso(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    #y_pred_Ridge = Ridge(ytrain = ytrain ,Xtrain = Xtrain, Xtest = Xtest)
    #y_pred_Ridge = MLP(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    
    # Evaluation
    corr_Ridge = cor(ytest,y_pred_Ridge)
    #print(corr_Ridge)
    #corr_RF = cor(ytest,y_pred_RF)
    #corr_ENet = cor(ytest,y_pred_ENet)
    #corr_Lasso = cor(ytest,y_pred_Lasso)
    #corr_Ridge = cor(ytest , y_pred_Ridge)
    #corr_Ridge = cor(ytest , y_pred_Ridge)
    
    ttest = t.test(y_pred_Ridge[ytest==1], y_pred_Ridge[ytest==2], alternative="greater")$p.value
    Ranksum = wilcox.test(y_pred_Ridge[ytest==1], y_pred_Ridge[ytest==2], alternative ="greater")$p.value
    #print(Ranksum)
  } else{
    corr_Ridge = 0
    ttest = 1
    Ranksum = 1
  }
    #plot(ytest,y_pred_Ridge)
    #boxplot(y_pred_Ridge[ytest==1], y_pred_Ridge[ytest==2])
    #plot(ytest,y_pred_SGL)
    #corr_MLP = cor(ytest,y_pred_MLP)
    result = data.frame(corr_Ridge = corr_Ridge, ttest=ttest, Ranksum = Ranksum)
                        #corr_SGL = corr_SGL,
                        #corr_RF = corr_RF,
                        #corr_ENet = corr_ENet,
                        #corr_Lasso = corr_Lasso,
                        #corr_MLP = corr_MLP)
    Results = rbind(Results, result)
    print(Results)
  }
    
sum(Results$Ranksum<0.05)
which(Results$Ranksum<0.05)
which(Results$ttest<0.05)
N_genes
# PRISM_TCGA_drugs = colnames(sen_PRISM)
# saveRDS(PRISM_TCGA_drugs,"Processed_data/Other/PRISM_TCGA_drugs.rds")
# which(!is.na(res_TCGA[,6]))
# boxplot(y_pred_Ridge[ytest==1], y_pred_Ridge[ytest==2],
#         y_pred_Ridge[ytest==3], y_pred_Ridge[ytest==4])
# 



