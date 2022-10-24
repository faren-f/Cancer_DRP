rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

GE_PRISM = readRDS("Processed_data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")
GE_TCGA = readRDS("Processed_data/S23/expresion_matrix_TCGA.rds")

# Remove genes whose Q3 is zero
q3_genes = apply(GE_TCGA,2,quantile,prob=0.75)
GE_TCGA = GE_TCGA[,-which(q3_genes==0)]
GE_PRISM = GE_PRISM[,-which(q3_genes==0)]

#write.csv(GE_PRISM,"~/Desktop/Cancer_DRP/Python/decouple_R/Raw_data/Data_from_R/GE_PRISM.csv")
#write.csv(GE_TCGA,"~/Desktop/Cancer_DRP/Python/decouple_R/Raw_data/Data_from_R/GE_TCGA.csv")
  
source("F15-Feature_Selection_PRISM@TCGA.R")
selected_features = c("TF_decoupleR")
Omics_List = Feature_Selection_PRISM_TCGA(selected_features, 
                                          Xtrain=GE_PRISM, 
                                          Xtest=GE_TCGA)
Xtrain = Omics_List[[1]]
index = Omics_List[[2]]
Xtest = Omics_List[[3]]

saveRDS(Omics_List, "Processed_data/S34/fgsea.rds")

