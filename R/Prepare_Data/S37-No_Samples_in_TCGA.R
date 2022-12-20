rm(list=ls())
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

res_TCGA = readRDS("Processed_data/S21/Drug_Response_matrix_TCGA.rds")

for(i in 1:ncol(res_TCGA)){
  res_TCGA[which(res_TCGA[,i]== 2),i]=1
  res_TCGA[which(res_TCGA[,i]== 3 | res_TCGA[,i]== 4),i]= 2
}

S1 = c()
S2 = c()
for(i in 1:ncol(res_TCGA)){
  S1 = c(S1, sum(res_TCGA[,i]== 1, na.rm = TRUE))
  S2 = c(S2, sum(res_TCGA[,i]== 2, na.rm = TRUE))
}
S1 = data.frame(S1)  
S2 = data.frame(S2)  
S = cbind(S1,S2)
rownames(S) = colnames(res_TCGA)

a = rownames(S)[which(S[,1]>5 & S[,2]>5)]

sen = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")
sen = sen[,-c(3,10,12,14)]
b = colnames(sen)

intersect(a,b)
