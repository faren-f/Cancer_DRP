rm(list=ls())

library(ggplot2)
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
Tissue = readRDS("Processed_data/S19/sample_tissue_types.rds")
Tissue_TCGA = readRDS("Processed_data/S21/Cancer_types_TCGA.rds")
N_Tissue_TCGA = readRDS("Processed_data/S22/Number_of_each_Cancer_TCGA.rds")
sen_PRISM = readRDS("Processed_data/S1/sensitivity_matrix_AUC.rds")
res_TCGA = readRDS("Processed_data/S23/Drug_response_matrix_TCGA.rds")

GE = readRDS("Processed_data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")
GE_TCGA = readRDS("Processed_data/S23/expresion_matrix_TCGA.rds")

GE_TCGA[is.na(GE_TCGA)] = 0
r = apply(GE_TCGA,2,quantile,prob=0.75)
GE_TCGA = GE_TCGA[,-which(r==0)]
GE = GE[,-which(r==0)]


# Sample_Tissue_matrix for TCGA
Sample_Tissue_mat = matrix(0,nrow(GE_TCGA),nrow(N_Tissue_TCGA))
rownames(Sample_Tissue_mat) = rownames(GE_TCGA)
colnames(Sample_Tissue_mat) = N_Tissue_TCGA[,1]

k=1
N = as.numeric(N_Tissue_TCGA[,2])
for(j in 1:ncol(Sample_Tissue_mat)){
  Sample_Tissue_mat[c(k:(N[j]+(k-1))),j]=1
  k = k+N[j]
} 

#I_G = intersect(l1000_genes,colnames(GE))
#GE = GE[,I_G]

#plot data
pca_PRISM = prcomp(GE, scale. = T)
pc_PRISM = pca_PRISM$x
plot(pc_PRISM[,1],pc_PRISM[,2])

pca_TCGA = prcomp(GE_TCGA, scale. = T)
pc_TCGA = pca_TCGA$x
plot(pc_TCGA[,1],pc_TCGA[,2])

# outlier detection
which(pc_TCGA[,1]< -500)
GE_TCGA = GE_TCGA[c(-194,-719),]
GE = GE[c(-194,-719),]
Sample_Tissue_mat = Sample_Tissue_mat[c(-194,-719),]
  
pca_TCGA = prcomp(GE_TCGA, scale. = T)
pc_TCGA = pca_TCGA$x
plot(pc_TCGA[,1],pc_TCGA[,2])


# plot gene expresion of cell lines with different colors for different tissues
j=0
T_i = matrix(0,nrow(GE),1)
rownames(T_i) = rownames(GE)
for (i in colnames(Tissue)){
  j=j+1
  p = which(Tissue[,i]==1)
  T_i[p] = j
}
Tissues = data.frame(colnames(Tissue))
# myCol = c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(7, "Accent"))
# plot(pc[,2],pc[,3], col = myCol[T_i],pch = 20)
plot(pc[,2],pc[,3], col = ifelse(Tissue[,4]==1, alpha("red",.8), alpha("grey",.3)), pch = 20)

# plot gene expresion of Tissues with different colors for different tissues

x = c(1:26)
x = data.frame(x,Rep = N_Tissue_TCGA[,2])
T_i_TCGA = data.frame(lapply(x, rep, x$Rep))
T_i_TCGA = T_i_TCGA[,1]
myCol = c(RColorBrewer::brewer.pal(12, "Set3"), 
          RColorBrewer::brewer.pal(11, "Paired"),
          RColorBrewer::brewer.pal(3, "Set1"))
plot(pc_TCGA[,2],pc_TCGA[,3], col = myCol[T_i_TCGA],pch = 20)

plot(pc_TCGA[,1],pc_TCGA[,2], col = ifelse(Sample_Tissue_mat[,17]==1, 
                                 alpha("red",1), alpha("grey",.3)), pch = 20)
Tissue_TCGA = data.frame(Tissue_TCGA)

# plot gene expresion of cell lines with different colors for different tissues 
#and diffrent size base on their sensitivies
sen = sen_PRISM[,1362]
sen[is.na(sen)] = 0
sen = sen/max(sen)
plot(pc[,2],pc[,3], cex = 2*sen,
     col = ifelse(Tissue[,21]==1, alpha("red",.8), alpha("grey",.3)), pch = 20)




#general codes for PCA analysis
#Q = pca$rotation
#biplot(pca, scale = 0)
#calculate total variance explained by each principal component
#var_explained = pca$sdev^2 / sum(pca$sdev^2)
#create scree plot
# qplot(c(1:ncol(GE)), var_explained) + 
#   geom_line() + 
#   xlab("Principal Component") + 
#   ylab("Variance Explained") +
#   ggtitle("Scree Plot")+
#   ylim(0,1)
