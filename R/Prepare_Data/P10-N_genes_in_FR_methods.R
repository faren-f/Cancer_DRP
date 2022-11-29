rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
sen_PRISM = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")

GE_PRISM = readRDS("Processed_data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")
GE_TCGA = readRDS("Processed_data/S23/expresion_matrix_TCGA.rds")

# Remove genes whose Q3 is zero
q3_genes = apply(GE_TCGA,2,quantile,prob=0.75)
GE_TCGA = GE_TCGA[,-which(q3_genes==0)]
GE_PRISM = GE_PRISM[,-which(q3_genes==0)]

#Landmark
l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")

#TF activity
dR_PRISM = read.table("Processed_data/S33/gsea2_PRISM.csv",sep = ",",header = TRUE, row.names = 1)

# Progeny
source("F15-Feature_Selection_PRISM@TCGA.R")
selected_features = c("progeny")
Omics_List = Feature_Selection_PRISM_TCGA(selected_features, Xtrain=GE_PRISM, Xtest=GE_TCGA)
GE_progeny = Omics_List[[1]]

#PW
source("F22-Drug_Pathway_Level_genes_eachTarget.R")
N_genes = c()
for(i in 1:ncol(sen_PRISM)){
  print(i)
  pathway_gene_set = Drug_Pathway_gene_set_eachTarget(drug = i, level=1)[["all"]]
  I = intersect(colnames(GE_PRISM),pathway_gene_set)
  N_genes = c(N_genes, length(I))
}


pdf(paste0("Figures/FS/Number_of_genes_in_FR_Methods/BoxPlot_N_genes_FR_Methods.pdf"),
    height = 4, width = 4)
boxplot(ncol(GE_PRISM), length(l1000_genes), N_genes[c(-22,-23)],names=NA,
        ncol(GE_progeny),
        ncol(dR_PRISM))
        #medcol=c("#F2E6D6","#A0E0E4","#FCC0C5","#A49393","#F582A8"),
        #col = c("#F2E6D6","#A0E0E4","#FCC0C5","#A49393","#F582A8"),
        #border = 4, boxcol = 2,   # Box border color
        #boxfill = "#FCC0C5")  # Box fill color)
dev.off()




