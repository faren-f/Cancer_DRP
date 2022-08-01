
#                      Created on Sun Jul 31 2022

#                     @author: Farzaneh Firoozbakht

# Discription: 
# This script reads ppi_EdgeList_compelete that is obtained from Step7, 
# gene expresion matrix that is obtained from Step1,
# to find the intersect gene symbols between ppi_EdgeList and GE.


rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

ppi_edgelist = readRDS("ppi_data/Processed_data/ppi_EdgeList_compelete.rds")
GE = readRDS("Processed_Data/expresion_matrix.rds")

#Removing the gene symbol rows with "NA"
ppi_edgelist = ppi_edgelist[!is.na(ppi_edgelist$gene_symbol1),]
ppi_edgelist = ppi_edgelist[!is.na(ppi_edgelist$gene_symbol2),]


# Since length(intersect(ppi_edgelist$gene_symbol1,ppi_edgelist$gene_symbol2)) =
#length(unique(ppi_edgelist$gene_symbol1)) = length(unique(ppi_edgelist$gene_symbol2)),
# intersect_GE_ppi is obtained based on just (gene_symbol1)

intersect_GE_ppi = intersect(colnames(GE),ppi_edgelist$gene_symbol1)
GE = GE[,intersect_GE_ppi]
ppi_edgelist = ppi_edgelist[ppi_edgelist$gene_symbol1 %in% intersect_GE_ppi,]
ppi_edgelist = ppi_edgelist[ppi_edgelist$gene_symbol2 %in% intersect_GE_ppi,]

saveRDS(ppi_edgelist,"ppi_data/Processed_data/ppi_EdgeList_compelete_PRISM.rds")
saveRDS(GE,"ppi_data/Processed_data/expresion_matrix_PRISM_ppi.rds")

