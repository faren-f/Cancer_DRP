rm(list=ls())
library("readxl")

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
TCGA_Table = read_excel("Raw_data/TCGA/2016-Ding/bioinfo16_supplementary_tables.xlsx",
                        sheet = 3,na = "---")

colnames(TCGA_Table) = TCGA_Table[2,]
TCGA_Table = TCGA_Table[c(-1,-2),]
cancers = unique(TCGA_Table$Cancer)
cancers = cancers[-1]

P = c()
for (i in cancers){
  P_i = TCGA_Table[TCGA_Table$Cancer==i,2]
  P_i = P_i[!is.na(P_i)]
  P = c(P,length(unique(P_i)))
  
}




drugs = unique(TCGA_Table$drug_name)
drugs = drugs[-1]
drugs = tolower(drugs)

sen = readRDS("Processed_data/S1/sensitivity_matrix_AUC.rds")
TCGA_PRISM_drugs = intersect(drugs,colnames(sen))
saveRDS(TCGA_PRISM_drugs,"Processed_data/S21/TCGA_PRISM_drugs.rds")

BLCA = read.table("Raw_data/TCGA/RNAseq/BLCA/BLCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",
                  fill = TRUE, header=FALSE)

S = t(data.frame(BLCA[1,]))
  
  
  


