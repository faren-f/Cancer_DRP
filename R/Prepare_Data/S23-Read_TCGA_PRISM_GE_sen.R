rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
drugs = readRDS("Processed_data/S21/Drugs_TCGA@PRISM.rds")
res_TCGA = readRDS("Processed_data/S22/Drug_response_matrix_TCGA.rds")
sen = readRDS("Processed_data/S1/sensitivity_matrix_AUC.rds")
sen_G = readRDS("All_Results/sen_PRISM_good_drugs.rds")

GE_PRISM = readRDS("Processed_data/S1/expresion_matrix.rds")
GE_TCGA = readRDS("Processed_data/S22/expresion_matrix_TCGA.rds")
I_GE = intersect(colnames(GE_PRISM),colnames(GE_TCGA))
GE_PRISM = GE_PRISM[,I_GE]
GE_TCGA = GE_TCGA[,I_GE]

# Intercept of good drugs in PRISM with TCGA
I_D_G = intersect(colnames(sen_G),drugs)
sen_PRISM_G = sen_G[,I_D_G]
res_TCGA_G = res_TCGA[,I_D_G]

# Intercept of all drugs in PRISM with TCGA
I_D = intersect(colnames(sen),drugs)
sen_PRISM = sen[,I_D]
res_TCGA = res_TCGA[,I_D]


# Save data
saveRDS(GE_PRISM,"Processed_data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")
saveRDS(GE_TCGA,"Processed_data/S23/expresion_matrix_TCGA.rds")

#saveRDS(sen_PRISM_G,"Processed_data/S23/sensitivity_matrix_PRISM_with@TCGA@GoodDrugs.rds")
#saveRDS(res_TCGA_G,"Processed_data/S23/Drug_response_matrix_TCGA@goodDrugs.rds")

saveRDS(sen_PRISM,"Processed_data/S23/sensitivity_matrix_PRISM_with@TCGA@drugs.rds")
saveRDS(res_TCGA,"Processed_data/S23/Drug_response_matrix_TCGA.rds")

