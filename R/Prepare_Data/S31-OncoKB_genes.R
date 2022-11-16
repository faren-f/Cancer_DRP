#                    Created on Sat Oct. 22 13:31 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This script receives OncoKB gene list from https://www.oncokb.org
# and find 1066 unique OncoKB genes

rm(list=ls())
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

OncoKB = read.table("Raw_data/OncoKB/cancerGeneList.tsv", sep = '\t', header = TRUE)

OncoKB_genes = OncoKB$Hugo.Symbol
OncoKB_oncogenes = OncoKB$Hugo.Symbol[OncoKB$Is.Oncogene=="Yes"]

saveRDS(OncoKB_genes, "Processed_data/S31/OncoKB_genes.rds")
saveRDS(OncoKB_oncogenes, "Processed_data/S31/OncoKB_oncogenes.rds")



