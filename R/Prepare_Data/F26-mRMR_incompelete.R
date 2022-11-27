rm(list=ls())

library(praznik)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
sen_PRISM = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")
res_TCGA = readRDS("Processed_data/Other/Res_TCGA_24_Drugs.rds")

GE_PRISM = readRDS("Processed_data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")
GE_TCGA = readRDS("Processed_data/S23/expresion_matrix_TCGA.rds")

# Remove genes whose Q3 is zero
q3_genes = apply(GE_TCGA,2,quantile,prob=0.75)
GE_TCGA = GE_TCGA[,-which(q3_genes==0)]
GE_PRISM = GE_PRISM[,-which(q3_genes==0)]

i=5
Xtrain = GE_PRISM[!is.na(sen_PRISM[,i]),]
ytrain = sen_PRISM[!is.na(sen_PRISM[,i]),i]

Xtrain = data.frame(Xtrain)

mRMR = MRMR(Xtrain, ytrain, k = 500)
plot(mRMR$score[1:10])


D = c()
for(j in 1:(length(mRMR$score)-2)){
  D1 = mRMR$score[j]-mRMR$score[j+1]
  D2 = mRMR$score[j+1]-mRMR$score[j+2]
  D = c(D, D2/D1)
}

plot(D)






