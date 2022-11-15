rm(list=ls())

library(parallel)
no_cores = detectCores()
cl = makeCluster(no_cores-2)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

sen_GDSC = readRDS("Processed_Data/S39/GDSC_sensitivity_matrix_AUC.rds")
colnames(sen_GDSC) = tolower(colnames(sen_GDSC))
res_TCGA = readRDS("Processed_data/S24/Drug_response_TCGA_binarized.rds")

intersected_drugs = intersect(colnames(sen_GDSC),colnames(res_TCGA))
sen_GDSC = sen_GDSC[,intersected_drugs]
res_TCGA = res_TCGA[,intersected_drugs]

GE_GDSC = readRDS("Processed_Data/S39/GDSC_expresion_matrix.rds")
GE_TCGA = readRDS("Processed_data/S23/expresion_matrix_TCGA.rds")

intersected_genes = intersect(colnames(GE_GDSC),colnames(GE_TCGA))
GE_GDSC = GE_GDSC[,intersected_genes]
GE_TCGA = GE_TCGA[,intersected_genes]

# Remove genes whose Q3 is zero
q3_genes = apply(GE_TCGA,2,quantile,prob=0.75)
GE_TCGA = GE_TCGA[,-which(q3_genes==0)]
GE_GDSC = GE_GDSC[,-which(q3_genes==0)]

clusterExport(cl, c("GE_GDSC","GE_TCGA","sen_GDSC","res_TCGA"))
clusterEvalQ(cl, c(source("F10-Ridge.R"), source("F18-Combat_Normalization.R")))

DrugLoop = function(i){
  
  print(paste0("The drug number is: ", as.character(i)))
  
  Xtrain = GE_GDSC[!is.na(sen_GDSC[,i]),]
  ytrain = sen_GDSC[!is.na(sen_GDSC[,i]),i]
  
  Xtest = GE_TCGA[!is.na(res_TCGA[,i]),]
  ytest = res_TCGA[!is.na(res_TCGA[,i]),i]
  
  if(sum(ytest==1)>3 & sum(ytest==2)>3){
    X_Normalization = Combat_Scale(Xtrain,Xtest)
    
    Xtrain = X_Normalization[[1]]
    Xtest = X_Normalization[[2]]
    
    source("F15-Feature_Selection_PRISM@TCGA.R")
    selected_features = c("Landmark_genes")
    Omics_List = Feature_Selection_PRISM_TCGA(selected_features, Xtrain=Xtrain, Xtest=Xtest)
    Xtrain = Omics_List[[1]]
    index = Omics_List[[2]]
    Xtest = Omics_List[[3]]
    
    # Ytrain normalization
    ytrain = scale(ytrain)
    ytrain = ytrain[,1]
    # Models
    y_pred = Ridge(ytrain = ytrain ,Xtrain = Xtrain, Xtest = Xtest)
    
    # Evaluation
    corr = cor(ytest,y_pred)
    ttest = t.test(y_pred[ytest==1], y_pred[ytest==2], alternative="greater")$p.value
    Ranksum = wilcox.test(y_pred[ytest==1], y_pred[ytest==2], alternative ="greater")$p.value
    result = data.frame(corr = corr, ttest=ttest, Ranksum = Ranksum)
    
    }else{
      corr = 0
      ttest = 1
      Ranksum = 1
      result = data.frame(corr = corr, ttest=ttest, Ranksum = Ranksum)
    }
  return(result)
}

N_drug = ncol(sen_GDSC)
result = parLapply(cl, sapply(1:N_drug, list), DrugLoop) 

Result = data.frame()
for (k in 1:N_drug){
  Result = rbind(Result, result[[k]])
}

stopCluster(cl)

saveRDS(Result,"Final_Result/Train@PRISM_Test@TCGA_FS/RF2-Landmark.rds")
print(sum(Result$Ranksum<0.05))
print(which(Result$Ranksum<0.05))
print(which(Result$ttest<0.05))

#boxplot(y_pred[ytest==1],y_pred[ytest==2])


