rm(list=ls())
source("F18-Combat_Normalization.R")
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

GE = readRDS("Processed_data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")
GE_TCGA = readRDS("Processed_data/S23/expresion_matrix_TCGA.rds")

# Remove genes whose Q3 is zero
q3 = apply(GE_TCGA,2,quantile,prob=0.75)
GE_TCGA = GE_TCGA[,-which(q3==0)]
GE = GE[,-which(q3==0)]

batches = rep(c('cell', 'clin'), times = c(nrow(GE), nrow(GE_TCGA)))
batches = as.factor(batches)

X_Normalization = Combat_Scale(Xtrain = GE,Xtest = GE_TCGA)
GE_after = X_Normalization[[1]]
GE_TCGA_after = X_Normalization[[2]]


#### PCA
n_pc = 2
PCA = prcomp(rbind(GE, GE_TCGA))
PC = as.matrix(PCA$x)
PC_before = as.matrix(PC[,1:n_pc])

PCA = prcomp(rbind(GE_after, GE_TCGA_after))
PC = as.matrix(PCA$x)
PC_after = as.matrix(PC[,1:n_pc])


plot(PC_before[,1], PC_before[,2], col = ifelse(batches=="clin","red","blue"),
     xlab="PC1", ylab="PC2", pch = 20, cex = .5)

plot(PC_after[,1], PC_after[,2], col = ifelse(batches=="clin","red","blue"),
     xlab="PC1", ylab="PC2", pch = 20, cex = .5)

## Save data

saveRDS(GE_after,"Processed_data/S25/GE_PRISM_AfterCombat_withTCGA.rds")
saveRDS(GE_TCGA_after,"Processed_data/S25/GE_TCGA_AfterCombat_withPRISM.rds")


