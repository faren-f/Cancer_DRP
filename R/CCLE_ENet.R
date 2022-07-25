rm(list = ls())

library(GEOquery)
require(caTools)
library(keras)
library(corrplot)
library(gelnet)
library(RColorBrewer)

setwd("~/Desktop/R_Root/DRP/Paper_No2/")

#gseID = "GSE36139"
#data = getGEO(gseID)[[1]]
#saveRDS(data, paste0("data/",gseID,".rds"))

## Read data
Raw_data = readRDS("data/CCLE/getGEO/mRNA/GSE36139.rds")
sensitivity_data = read.csv("data/CCLE/CCLE_web/Pharmacological Profile/CCLE_NP24.2009_Drug_data_2015.02.24.csv")

featureData = fData(Raw_data)
phenoData = pData(Raw_data)
exprData_all = t(exprs(Raw_data))

## Sensitivity Matrix All
No.drugs = length(unique(sensitivity_data$Compound))

#sensitivity matrix for all the cell lines(some of them do not have gene expression)
sen_mat_all = matrix(0,length(unique(sensitivity_data$Primary.Cell.Line.Name)),No.drugs)

rownames(sen_mat_all) = unique(sensitivity_data$Primary.Cell.Line.Name)
colnames(sen_mat_all) = unique(sensitivity_data$Compound)

for (i in rownames(sen_mat_all)){
  for(j in colnames(sen_mat_all)){
    
    if (sum(sensitivity_data$Primary.Cell.Line.Name == i & sensitivity_data$Compound == j) == 0){
      sen_mat_all[i,j] = NA}
    else{
      sen_mat_all[i,j] = as.numeric(sensitivity_data$ActArea[sensitivity_data$Primary.Cell.Line.Name == i 
                                                             & sensitivity_data$Compound == j]) 
    }}
}
sen_mat_all = as.data.frame(sen_mat_all)
#plot(sort(sen_mat_all[,24]),pch = 1, cex = .3)

## Expression Matrix
row.names(exprData_all) = phenoData$title 
exprData = exprData_all[rownames(exprData_all) %in% rownames(sen_mat_all),]
No.celllines = nrow(exprData)

## Sensitivity Matrix [For Cell lines with Expression]
sen_mat = sen_mat_all[rownames(sen_mat_all) %in% rownames(exprData),]

# Match the rownames in sensitivity matrix & expression matrix 
sen_mat = sen_mat[rownames(exprData),]

## Pre Processing
# Remove Genes with Low Median
median_genes = apply(exprData,2,median)
names(median_genes) = colnames(exprData)
median_genes = sort(median_genes)
#hist(median_genes,100)
#abline(v=3.5,col="red")
median_genes = median_genes[median_genes>2]
exprData = exprData[,names(median_genes)]

# Remove Genes with Low STD
std_genes = apply(exprData,2,sd)
names(std_genes) = colnames(exprData)
std_genes = sort(std_genes)
hist(std_genes,100)
abline(v=.4,col="red")
std_genes = std_genes[std_genes>.3]
exprData = exprData[,names(std_genes)]

## Classifier
Final_Cor = rep(0,20)
mse_test = rep(0,20)

for (j in 1:20){
  i=5
  print(j)
  #for(i in 1:ncol(sen_mat)){
  sen_each_drug = sen_mat[,i]
  sen_each_drug = sen_each_drug[!is.na(sen_each_drug)]
  expr_each_drug = exprData[!is.na(sen_mat[,i]),] # For each drug some of the cell lines are "NA" so hear we remove them 
  
  ## Corr input & output
  Corr = cor(expr_each_drug,sen_each_drug)  
  #max(Corr)  
  hist(abs(Corr))  
  ind_Corr = which(abs(Corr)> 0.1)
  sel_expr_each_drug = expr_each_drug[,ind_Corr]
  #expr_each_norm = scale(sel_expr_each_drug)
  num = sel_expr_each_drug-min(sel_expr_each_drug)
  denum = max(sel_expr_each_drug)-min(sel_expr_each_drug)
  expr_each_drug_norm = num/denum
  
  #sen/res
  #median_sen_mat = median(sen_each_drug)
  #sen_each_drug_binarized = as.factor(ifelse (sen_each_drug>median_sen_mat,0,1))
  sen_each_drug_binarized = sen_each_drug
  num = sen_each_drug_binarized-min(sen_each_drug_binarized)
  denum = max(sen_each_drug_binarized)-min(sen_each_drug_binarized)
  sen_each_drug_binarized_norm = num/denum
  
  ## Split data into train & test
  sample = sample.split(sen_each_drug_binarized_norm, SplitRatio = .8)
  
  Xtrain = subset(expr_each_drug_norm, sample == TRUE)
  Xtest  = subset(expr_each_drug_norm, sample == FALSE)
  Ytrain = subset(sen_each_drug_binarized_norm, sample == TRUE)
  Ytest  = subset(sen_each_drug_binarized_norm, sample == FALSE)
  
  
  ## Model training
  n_feat = ncol(expr_each_drug_norm)
  lambda1 = 0.01
  lambda2 = 20
  d = rep(1, n_feat)
  
  GelNet = gelnet(Xtrain, Ytrain, l1 = lambda1, l2 = lambda2, d = d,
                  P = diag(d), m = rep(0,n_feat), max.iter = 50, eps = 1e-05)
  
  ## Test the model
  Beta = GelNet[["w"]]
  Beta0 = GelNet[["b"]]
  y_hat = (Xtest %*%  Beta) + Beta0
  y = Ytest
  
  mse_test[j] = mean((y - y_hat)^2)
  Final_Cor[j] = cor(y,y_hat, method = "pearson")

   print(mse_test)
   print(cor(y,y_hat))
}
#}

print(mean(mse_test))
print(mean(Final_Cor))
print(sd(Final_Cor))
