#                    Created on Sat Oct. 22 15:42  2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This script receives Cancer genes(CGs) from the supplimentary table S2A 
#(Iorio et al, 2016) and find 471 unique genes.

rm(list=ls())
library("readxl")

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
Cancer_genes_GDSC = read_excel("Raw_data/Cancer_genes_GDSC/1-s2.0-S0092867416307462-mmc3.xlsx",
                            sheet = 1)
colnames(Cancer_genes_GDSC) = Cancer_genes_GDSC[1,]
Cancer_genes_GDSC = Cancer_genes_GDSC[-1,]
Cancer_genes_GDSC = data.frame(Cancer_genes_GDSC)
Cancer_genes_GDSC = unique(Cancer_genes_GDSC[,1])
Cancer_genes_GDSC = Cancer_genes_GDSC[1:470]
saveRDS(Cancer_genes_GDSC, "Processed_data/S33/CancerGenesGDSC.rds")


