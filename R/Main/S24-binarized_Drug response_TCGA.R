rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
#res_TCGA_G = readRDS("Processed_data/S23/Drug_response_matrix_TCGA@goodDrugs.rds")
res_TCGA = readRDS("Processed_data/S23/Drug_response_matrix_TCGA.rds")

N_Cancer = readRDS("Processed_data/S22/Number_of_each_Cancer_TCGA.rds")

GE = readRDS("Processed_data/S23/expresion_matrix_TCGA.rds")
GE[is.na(GE)] = 0

N_Cancer = data.frame(N_Cancer)
N_Cancer[,2] = as.numeric(N_Cancer[,2])

T_i = data.frame(lapply(N_Cancer, rep, N_Cancer$N_TT))
GE = GE[T_i$N_TT>3,]

for(i in 1:ncol(res_TCGA)){
  res_TCGA[which(res_TCGA[,i]== 2),i]=1
  res_TCGA[which(res_TCGA[,i]== 3 | res_TCGA[,i]== 4),i]= 2
}

saveRDS(res_TCGA,"Processed_data/S24/Drug_response_TCGA_binarized.rds")
#saveRDS(res_TCGA_G,"Processed_data/S24/Drug_response_TCGA_binarized.rds")





