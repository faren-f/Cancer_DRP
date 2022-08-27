rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

GE_PRISM = readRDS("Processed_data/S1/expresion_matrix.rds")
GE_TCGA = readRDS("Processed_data/S22/expresion_matrix_TCGA.rds")
I_GE = intersect(colnames(GE_PRISM),colnames(GE_TCGA))
GE_PRISM = GE_PRISM[,I_GE]
GE_TCGA = GE_TCGA[,I_GE]

drugs = readRDS("Processed_data/S21/Drugs_TCGA@PRISM.rds")
sen = readRDS("All_Results/sen_PRISM_good_drugs.rds")
I_D = intersect(colnames(sen),drugs)
sen_PRISM = sen[,I_D]

res_TCGA = readRDS("Processed_data/S22/Drug_response_matrix_TCGA.rds")
res_TCGA = res_TCGA[,I_D]

# Save data
saveRDS(GE_PRISM,"Processed_data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")
saveRDS(GE_TCGA,"Processed_data/S23/expresion_matrix_TCGA.rds")

saveRDS(sen_PRISM,"Processed_data/S23/sensitivity_matrix_PRISM_with@TCGA@drugs.rds")
saveRDS(res_TCGA,"Processed_data/S23/Drug_response_matrix_TCGA.rds")

