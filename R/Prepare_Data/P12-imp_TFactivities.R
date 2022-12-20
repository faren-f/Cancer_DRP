rm(list=ls())

library(ROCR)
source("F18-Combat_Normalization.R")
library(glmnet)
library(caret)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

dR_PRISM = read.table("Processed_data/S33/gsea2_PRISM.csv",sep = ",",header = TRUE, row.names = 1)
dR_TCGA = read.table("Processed_data/S33/gsea2_TCGA.csv",sep = ",",header = TRUE, row.names = 1)
Cor = cor(dR_TCGA)
diag(Cor) = 0
Max_Cor_TFs = apply(Cor,2,max)
hist(Max_Cor_TFs)

sen_PRISM = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")
res_TCGA = readRDS("Processed_data/Other/Res_TCGA_24_Drugs.rds")

# Remove genes whose Q3 is zero
q3_genes = apply(dR_TCGA,2,quantile,prob=0.75)
if(sum(q3_genes==0)>0){
  dR_TCGA = dR_TCGA[,-which(q3_genes==0)]
  dR_PRISM = dR_PRISM[,-which(q3_genes==0)]
}

SigDrugs = c("etoposide","paclitaxel","leucovorin", 
             "ifosfamide", "gemcitabine",
             "cisplatin", "vinblastine")

result = c()
Wilcox_Test_AllDrugs = c()
T_Test_AllDrugs = c()

#for (i in SigDrugs){
  print(paste0("The drug number is: ", i))
  
  Xtrain = dR_PRISM[!is.na(sen_PRISM[,i]),]
  ytrain = sen_PRISM[!is.na(sen_PRISM[,i]),i]
  
  Xtest = dR_TCGA[!is.na(res_TCGA[,i]),]
  ytest = res_TCGA[!is.na(res_TCGA[,i]),i]
  
  X_Normalization = Combat_Scale(Xtrain,Xtest)
  
  Xtrain = X_Normalization[[1]]
  Xtest = X_Normalization[[2]]
  
  # Ytrain normalization
  ytrain = scale(ytrain)
  ytrain = ytrain[,1]
  
  train_data = cbind(Xtrain,ytrain)
  control = trainControl(method = "repeatedcv",
                         number = 5,
                         repeats = 5,
                         verboseIter = FALSE)
  
  #tune = expand.grid(alpha = 0,lambda = seq(0.01,5,by = 0.01))
  tune = expand.grid(alpha = 0,lambda = seq(0.01,5,by = 0.1))

  model = caret::train(ytrain ~., data = train_data,
                       method = "glmnet",
                       metric="RMSE",
                       allowParallel = TRUE,
                       tuneGrid = tune,
                       trControl = control)
  y_pred = predict(model,Xtest)
  Ranksum = wilcox.test(y_pred[ytest==1], y_pred[ytest==2], alternative ="greater")$p.value
  
  Beta = as.matrix(coef(model$finalModel, model$bestTune$lambda))

  Beta_Null = c()
  Wilcox_Test = c()
  T_test = c()
  # for(j in 1:1000){
  #   print(j)
  #   ytrain_perm = sample(ytrain)
  #   train_data = cbind(Xtrain, ytrain_perm)
  #   control = trainControl(method = "repeatedcv",
  #                          number = 5,
  #                          repeats = 5,
  #                          verboseIter = FALSE)
  #   
  #   #tune = expand.grid(alpha = 0,lambda = seq(0.01,5,by = 0.05))
  #   tune = expand.grid(alpha = 0, lambda = seq(0.01,5,by = 0.1))
  #   
  #   model = caret::train(ytrain_perm ~., data = train_data,
  #                        method = "glmnet",
  #                        metric="RMSE",
  #                        allowParallel = TRUE,
  #                        tuneGrid = tune,
  #                        trControl = control)
  #   y_pred = predict(model,Xtest)
  #   beta = as.matrix(coef(model$finalModel, model$bestTune$lambda))
  #   Beta_Null = cbind(Beta_Null, beta)
  # }
  
  
  #saveRDS(Beta_Null, "Processed_data/P12/Beta_Null_ifosfamide.rds")
  
  Beta_Null = readRDS("Processed_data/P12/Beta_Null_etoposide.rds")
  # Beta_Null2 = Beta_Null
  # Beta_Null = cbind(Beta_Null1,Beta_Null2)
  
  #P_val_Rank = c()
  PNorm = c()
  #p.val_KS = c()
  source("F27-Probability_Norm.R")
  for(k in 1:nrow(Beta_Null)){
    
    P = Probability_norm(abs(Beta[k]), abs(Beta_Null[k,]), alternative = "greater")
    PNorm = c(PNorm, P)
    
    # p.val_KS = c(p.val_KS, ks.test(Beta_Null[k,], "pnorm", mean(Beta_Null[k,]), 
    #                      sd(Beta_Null[k,]), exact = TRUE)$p.value)
    # 
    #P = rank.test(abs(Beta_Null[k,]), abs(Beta[k]), alternative = "greater")
    #P_val_Rank = c(P_val_Rank, P)
  }
  # hist(p.val_KS,100)
  # which(p.val_KS<0.05)
  # p.val_KS[p.val_KS<0.05]
  
  hist(PNorm,100)
  hist(p.adjust(PNorm, method = "BH"), 100, xlim = c(0,1))
  p_adjust_PNorm = p.adjust(PNorm, method = "BH")
  
  #hist(p.adjust(P_val_Rank, method = "BH"), 100, xlim = c(0,1))
  #p_adjust_Rank = p.adjust(P_val_Rank, method = "none")
  
  #Wilcox_Test_AllDrugs = cbind(Wilcox_Test_AllDrugs,Wilcox_Test)
#}

a = which(PNorm<0.01)
PNorm[PNorm<0.01]


sum(p_adjust_PNorm<0.01)
a = which(p_adjust_PNorm<0.01)
p_adjust_PNorm[p_adjust_PNorm<0.01]


a = order(p_adjust_PNorm, decreasing = FALSE)[1:10]
sort(p_adjust_PNorm, decreasing = FALSE)[1:10]



TF_sig = rownames(Beta)[a]


# sum(p_adjust_Rank<0.01)
# which(p_adjust_Rank<0.01)
# p_adjust_Rank[p_adjust_Rank<0.01]

hist(Beta_Null[1,],100)

