#                    Created on Wed Aug 11 11:14 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This script receives gene expression data and dorothea package to 
# find transcription factors

rm(list=ls())

library(dorothea)
data("dorothea_hs", package = "dorothea")

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
X = readRDS("Processed_Data/S1/expresion_matrix.rds")
X = t(X)

TF = dorothea::run_viper(X, dorothea_hs,
                           options =  list(minsize = 5,
                                           eset.filter = FALSE, cores = 1,
                                           verbose = FALSE))

saveRDS(TF_activities,"Processed_data/S14/TF.rds")


