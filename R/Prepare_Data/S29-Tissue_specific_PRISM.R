rm(list=ls())
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

sen = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")
TCGA_Patients = readRDS("Processed_data/S22/TCGA_Patients.rds")
N_Cancer_TCGA = readRDS("Processed_data/S22/Number_of_each_Cancer_TCGA.rds")

sample_tissue_one_hot = readRDS("Processed_data/S19/sample_tissue_types.rds")

Cellline_Tissue = matrix(0,nrow(sample_tissue_one_hot),2)
rownames(Cellline_Tissue) = rownames(sample_tissue_one_hot)

for(i in 1:ncol(sample_tissue_one_hot)){
  Cellline_Tissue[sample_tissue_one_hot[,i]==1,1] = i
  Cellline_Tissue[sample_tissue_one_hot[,i]==1,2] = colnames(sample_tissue_one_hot)[i]
}
Cellline_Tissue = Cellline_Tissue[rownames(sen),]

# Removing "fibroblast" and "rhabdomyosarcoma" samples from Cellline_Tissue and sensitivity matrix
Cellline_Tissue = Cellline_Tissue[-which(Cellline_Tissue[,1] == 23),] 
Cellline_Tissue = Cellline_Tissue[-which(Cellline_Tissue[,1] == 24),] 
sen = sen[rownames(Cellline_Tissue),] 

# Removing "fibroblast" and "rhabdomyosarcoma" tissues from sample_tissue_one_hot
sample_tissue_one_hot = sample_tissue_one_hot[,c(-23,-24)]

Tissue_Sen = matrix(0,ncol(sample_tissue_one_hot),ncol(sen))
rownames(Tissue_Sen) = colnames(sample_tissue_one_hot)
colnames(Tissue_Sen) = colnames(sen)

for(i in 1:ncol(sample_tissue_one_hot)){
  Tissue_Sen[i,] = apply(sen[Cellline_Tissue[,1] %in% i,],2,sum,na.rm = TRUE)
}


dists = as.dist(1-cor(t(Tissue_Sen)))
Tree = hclust(dists,method="ward.D")
plot(Tree)
k=4
clusters = cutree(Tree, k=k)
table(clusters)


col.pal = RColorBrewer::brewer.pal(9, "Reds")
pheatmap::pheatmap(dists,show_rownames=F, show_colnames=F, treeheight_row = 0, 
                   clustering_method = "ward.D", cutree_cols=k)


pheatmap::pheatmap(1-cor(t(Tissue_Sen))[Tree$order,Tree$order],
                   cluster_row = F,
                   cluster_cols = F,
                   show_rownames=F, show_colnames=F,
                   fontsize = 6.5, fontsize_row=6,fontsize_col = 6)

pheatmap::pheatmap(dists,show_rownames=FALSE, show_colnames=FALSE, treeheight_row = 0,
                   cutree_cols=k, cutree_rows=k, clustering_method = "ward.D")  
  
  