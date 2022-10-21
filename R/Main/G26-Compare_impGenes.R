rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

GE_PRISM = readRDS("Processed_data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")
GE_TCGA = readRDS("Processed_data/S23/expresion_matrix_TCGA.rds")


#Result_PRISM = readRDS("Final_Result/imp_genes_PRISM&TCGA/Docetaxel/Result_PRISM.rds")
Order_Beta_PRISM = readRDS("Final_Result/imp_genes_PRISM&TCGA/Docetaxel/Order_Beta_PRISM_WholeGenes.rds")


#Result_TCGA = readRDS("Final_Result/imp_genes_PRISM&TCGA/Docetaxel/Result_TCGA.rds")
Order_Beta_TCGA = readRDS("Final_Result/imp_genes_PRISM&TCGA/Docetaxel/Order_Beta_TCGA_WholeGenes.rds")
RankSum_Beta = apply(Order_Beta_TCGA,1,sum)
plot(sort(RankSum_Beta), pch=20)

len_intersect = c()
for(i in seq(100,3000,by = 100)){
  
  Order_PRISM = Order_Beta_PRISM[1:i]
  Order_TCGA = order(RankSum_Beta)[1:i]
  len_intersect = c(len_intersect, length(intersect(Order_PRISM, Order_TCGA)))
}


PRISM_goodGenes = colnames(GE_PRISM)[Order_Beta_PRISM[2:51]]
TCGA_goodGenes = colnames(GE_TCGA)[order(RankSum_Beta)[2:51]]

write.table(PRISM_goodGenes,"Final_Result/imp_genes_PRISM&TCGA/Docetaxel/PRISM_goodGenes@50.csv",
          quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(TCGA_goodGenes,"Final_Result/imp_genes_PRISM&TCGA/Docetaxel/TCGA_goodGenes@50.csv",
          quote = FALSE, row.names = FALSE, col.names = FALSE)

