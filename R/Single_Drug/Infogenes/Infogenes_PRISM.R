rm(list= ls())

library(gelnet)
library(igraph)
library(pso)
require(caTools)

setwd("~/Desktop/Cancer_DRP/R/Single_Drug/Infogenes/")

GE = readRDS("Raw_data/PRISM/expresion_matrix_ppi.rds")
sen = readRDS("Raw_data/PRISM/sensitivity_matrix.rds")
ppi_edgelist = readRDS("Raw_data/PRISM/ppi_EdgeList_compelete_PRISM.rds")



# Build graph -------------------------------------------------------------
ppi = rbind(ppi_edgelist$gene_symbol1, ppi_edgelist$gene_symbol2)
MyGraph = simplify(graph(ppi, directed = FALSE))



# Remove cell lines that do not have drug response from GE and sen
i = 1                            # drug number
X = GE[!is.na(sen[,i]),]
Y = sen[!is.na(sen[,i]),i]

# Normalization
X = scale(X)
Y = scale(Y)

Rep = 1
for (i in 1:Rep){
  
  ## Split data into train & test
  sample = sample.split(Y, SplitRatio = .8)
  
  Xtr_val = subset(X, sample == TRUE)
  Xtest  = subset(X, sample == FALSE)
  Ytr_val = subset(Y, sample == TRUE)
  Ytest  = subset(Y, sample == FALSE)
  
  sample = sample.split(Ytr_val, SplitRatio = .8)
  
  Xtrain = subset(Xtr_val, sample == TRUE)
  Xval  = subset(Xtr_val, sample == FALSE)
  Ytrain = subset(Ytr_val, sample == TRUE)
  Yval  = subset(Ytr_val, sample == FALSE)
  
  
  
  
  
}


  