
#                      Created on Sun Aug 9 2022

#                     @author: Farzaneh Firoozbakht

# Discription: 
# This script reads ppi_Omnipath that is obtained from Step8, 
# gene expresion matrix that is obtained from Step1,
# to find the intersect gene symbols between them.

rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

ppi_edgelist = readRDS("Processed_data/Step8/ppi_Omnipath.rds")
GE = readRDS("Processed_Data/Step1/expresion_matrix.rds")


ppi_edgelist = ppi_edgelist[(ppi_edgelist$source_genesymbol %in% colnames(GE) & 
                               ppi_edgelist$target_genesymbol %in% colnames(GE)),]
colnames(ppi_edgelist) = c("gene_symbol1", "gene_symbol2")

ppi_edgelist_dataframe = data.frame(ppi_edgelist)
ppi_edgelist_dataframe = c(ppi_edgelist_dataframe[,1],ppi_edgelist_dataframe[,2])

GE = GE[,intersect(colnames(GE),unique(ppi_edgelist_dataframe))]


saveRDS(ppi_edgelist,"Processed_Data/Step10/ppi_Omnipath_PRISM.rds")
saveRDS(GE,"Processed_data/Step10/expresion_matrix_PRISM_Omnipath.rds")

