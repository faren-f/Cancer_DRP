rm(list=ls())
library(dplyr)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
LINCS_genes = read.table("Raw_data/LINCS/GSE92742_Broad_LINCS_gene_info_delta_landmark.txt",
                         fill = TRUE , header = FALSE)
LINCS_genes = mutate_all(LINCS_genes,na_if,"")
f = data.frame(unique(LINCS_genes$V2))
colnames(f) = f[1,1]
f = f[!is.na(f$pr_gene_symbol),]
f=f[!grepl(",", f)]
f=f[!grepl(")", f)]
f=f[!grepl("^[a-z ]+$", f)]


Landmark_genes = f[c(-1,-69,-98,-382,-462,-480,-486,-569,-914,-938,-975,-983)]

ff= data.frame(Landmark_genes)
saveRDS(Landmark_genes,"Processed_data/S18/Landmark_genes.rds")

