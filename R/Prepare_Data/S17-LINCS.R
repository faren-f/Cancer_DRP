rm(list=ls())
library(dplyr)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
LINCS_genes = read.table("Raw_data/LINCS/GSE92742_Broad_LINCS_gene_info_delta_landmark.txt",
                         fill = TRUE , header = FALSE)
LINCS_genes = mutate_all(LINCS_genes,na_if,"")
f = data.frame(unique(LINCS_genes$V2))
f = na.omit(f)
colnames(f) = f[1,1]
f = f[-1,]

f=f[!grepl(",", f)]
f=f[!grepl(")", f)]
f=f[!grepl("^[a-z ]+$", f)]
ff = data.frame(f)
ff = f[c(-68,-97,-381,-461,-479,-485,-568,-913,-937,-974,-982)]
saveRDS("Processed_data/")

