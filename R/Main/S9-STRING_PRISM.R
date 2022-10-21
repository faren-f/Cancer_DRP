
#                      Created on Sun Jul 31 2022

#                     @author: Farzaneh Firoozbakht

# Discription: 
# This script reads ppi_STRING that is obtained from Step7, 
# gene expresion matrix that is obtained from Step1,
# to find the intersect gene symbols between them.

rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

ppi_edgelist = readRDS("Processed_data/S7/ppi_STRING.rds")
GE = readRDS("Processed_Data/S1/expresion_matrix.rds")

#Removing the gene symbol rows with "NA"
ppi_edgelist = ppi_edgelist[!is.na(ppi_edgelist$gene_symbol1),]
ppi_edgelist = ppi_edgelist[!is.na(ppi_edgelist$gene_symbol2),]


ppi_edgelist = ppi_edgelist[(ppi_edgelist$gene_symbol1 %in% colnames(GE) & 
                               ppi_edgelist$gene_symbol2 %in% colnames(GE)),]
GE = GE[,intersect(colnames(GE),unique(c(ppi_edgelist[,5],ppi_edgelist[,6])))]


saveRDS(ppi_edgelist,"Processed_Data/S9/ppi_STRING_PRISM.rds")
saveRDS(GE,"Processed_data/S9/expresion_matrix_PRISM_STRING.rds")







