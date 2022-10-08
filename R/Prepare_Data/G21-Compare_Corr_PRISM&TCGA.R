rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
source("F14-Feature_Selection.R")
source("F15-Feature_Selection_PRISM@TCGA.R")
source("F18-Combat_Normalization.R")

sen_PRISM = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")
res_TCGA = readRDS("Processed_data/Other/Res_TCGA_24_Drugs.rds")

GE_PRISM = readRDS("Processed_data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")
GE_TCGA = readRDS("Processed_data/S23/expresion_matrix_TCGA.rds")

# Remove genes whose Q3 is zero
q3_genes = apply(GE_TCGA,2,quantile,prob=0.75)
GE_TCGA = GE_TCGA[,-which(q3_genes==0)]
GE_PRISM = GE_PRISM[,-which(q3_genes==0)]

inteval = seq(200,800,by = 200)
Intersect = matrix(0,ncol(sen_PRISM),length(inteval))
for (i in 1:ncol(sen_PRISM)){
  X_PRISM = GE_PRISM[!is.na(sen_PRISM[,i]),]
  y_PRISM = sen_PRISM[!is.na(sen_PRISM[,i]),i]
  
  X_TCGA = GE_TCGA[!is.na(res_TCGA[,i]),]
  y_TCGA = res_TCGA[!is.na(res_TCGA[,i]),i]
  
  X_Normalization = Combat_Scale(X_PRISM,X_TCGA)
  
  X_PRISM = X_Normalization[[1]]
  X_TCGA = X_Normalization[[2]]
  
  selected_features = c("Landmark_genes")
  Omics_List = Feature_Selection_PRISM_TCGA(selected_features, 
                                            Xtrain=X_PRISM, Xtest=X_TCGA)
  X_PRISM = Omics_List[[1]]
  X_TCGA = Omics_List[[3]]
  Corr_PRISM = cor(X_PRISM,y_PRISM)
  sort_Corr_PRISM = sort(Corr_PRISM,decreasing = TRUE)[1:20]
  Corr_TCGA = cor(X_TCGA,y_TCGA)
  sort_Corr_TCGA = sort(Corr_TCGA,decreasing = TRUE)[1:20]
  
  Intersect_i = c()
  for(j in inteval){
    order_Corr_PRISM = order(Corr_PRISM,decreasing = TRUE)[1:j]
    order_Corr_TCGA = order(Corr_TCGA,decreasing = TRUE)[1:j]
    Intersect_i = c(Intersect_i, length(intersect(order_Corr_PRISM, order_Corr_TCGA)))
  }
  Intersect[i,] = Intersect_i
}
rownames(Intersect) = colnames(sen_PRISM)
