rm(list= ls())

library(gelnet)
library(igraph)
library(pso)

setwd("~/Desktop/Cancer_DRP/R/Single_Drug/Infogenes/")

GE = readRDS("Raw_data/PRISM/expresion_matrix_ppi.rds")
sen = readRDS("Raw_data/PRISM/sensitivity_matrix.rds")
ppi_edgelist = readRDS("Raw_data/PRISM/ppi_EdgeList_compelete_PRISM.rds")



# Build graph -------------------------------------------------------------
ppi = rbind(ppi_edgelist$gene_symbol1, ppi_edgelist$gene_symbol2)
MyGraph = simplify(graph(ppi, directed = FALSE))






GE = GE[!is.na(sen[,i]),]
sen_i = sen[!is.na(sen[,i]),i]



  