rm(list=ls())

library(reactome.db)
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

conv_table = readRDS("Processed_data/S7/biomart_conversion_table.rds")
#Drugs = as.character(readRDS("Processed_data/Other/24_drugs.rds")[,1])
Drugs = colnames(readRDS("Processed_data/S1/sensitivity_matrix_AUC.rds"))

Drug_Targets = readRDS("Processed_data/S1/drug_targets.rds")
rownames(Drug_Targets) = Drug_Targets$name  

gene2path = as.list(reactomeEXTID2PATHID)
#path2gene = as.list(reactomePATHID2EXTID)
  
Pathways = list()
d = 0
for(i in Drugs){
  print(d)
  d= d+1
  
  Pathways[[i]] = list() 
  Drug_Targets_i = strsplit(Drug_Targets[i,2],", ")
  
  for(j in Drug_Targets_i[[1]]){
    entrez = conv_table[which(conv_table$hgnc_symbol%in% j)[1],4]
    pw = gene2path[[as.character(entrez)]]
    Pathways[[i]][[j]] = pw
  }
}

# Save data
saveRDS(Pathways,"Processed_data/S26/All_Drug_Pathways.rds")


