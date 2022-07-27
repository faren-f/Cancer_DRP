rm(list = ls())
setwd("~/Desktop/Codes/Cancer_DRP/R/")

# Library -----------------------------------------------------------------
library('rtracklayer')
library(ggplot2)
# Read Data ---------------------------------------------------------------
cellline_info = read.csv("Data/PRISM_Raw_Dataset/Secondary/secondary-screen-cell-line-info.csv")

response = read.csv("Data/PRISM_Raw_Dataset/Secondary/secondary-screen-dose-response-curve-parameters.csv")

RNAseq = read.table("Data/PRISM_Raw_Dataset/Expression/RNA_seq/CCLE_RNAseq_rsem_transcripts_tpm_20180929.txt.gz",
                    header = TRUE, check.names = FALSE)

#RNAseq = read.table("Data/PRISM_Raw_Dataset/Expression/RNA_seq/CCLE_RNAseq_genes_rpkm_20180929.gct.txt", 
                     #skip = 2, header = TRUE, sep = "\t")
                  

gene_transfer = import("Data/PRISM_Raw_Dataset/Expression/RNA_seq/gencode.v19.genes.v7_model.patched_contigs.gtf.gz")

gene_transfer = data.frame(gene_transfer)

# Pre-processing ----------------------------------------------------------
## Log normalization of genes
expr_raw = RNAseq[,c(-1,-2)]
expr = log2(expr_raw + 1)

### Depict Figure 1.A
#expr_mean = apply(expr, 1, mean)
#pdf(file = "Figures/Fig1_A_expr_mean.pdf", width = 5, height = 3.5)
#hist(expr_mean, 100,main = "Histogram of expression mean",xlab = " ")
#dev.off()

### Depict Figure 1.B
#expr_sd = apply(expr, 1, sd)
#pdf(file = "Figures/Fig1_B_expr_sd.pdf", width = 5, height = 3.5)
#hist(expr_sd, 100, xlim = c(0,2), main = "Histogram of expression SD", xlab = " ")
#dev.off()

expr = cbind(RNAseq[,1],expr)

## Remove expressions with low mean 
mean_expr = apply(expr[,-1], 1, mean)
expr = expr[mean_expr > 1,]

### Depict Figure 2.A
#mean_expr = mean_expr[mean_expr>1]
#pdf(file = "Figures/Fig2_A_expr_mean_filter.pdf", width = 5, height = 3.5)
#hist(mean_expr, 100, main = "Histogram of expression mean",xlab = " ")
#dev.off()

## Remove expressions with low std 
sd_expr = apply(expr[,-1], 1, sd)
expr = expr[sd_expr > 0.5,]

### Depict Figure 2.B
#sd_expr = sd_expr[sd_expr>0.5]
#pdf(file = "Figures/Fig2_B_expr_sd_filter.pdf", width = 5, height = 3.5)
#hist(sd_expr, 100, xlim = c(0,2), main = "Histogram of expression SD",xlab = " ")
#dev.off()


## Some of the numbers has been removed in the previous steps thus we 
                #assign new numbers to rownames

rownames(expr) = 1:nrow(expr)

## Remove cell lines that do not exist in response from expression
expr1 = expr[,-1]
ccle_name_intersect = colnames(expr1) %in% response$ccle_name
expr1 = expr1[,ccle_name_intersect]

## Convert ccle_name to depmap_id in expression matrix
id_name = response[,c(2,3)]
id_name = id_name[!duplicated(id_name[,2]),]
id_name = id_name[!is.na(id_name[,2]),]
rownames(id_name) = id_name[,2]
dep_id = id_name[colnames(expr1),1] 
colnames(expr1) = dep_id
expr = cbind(expr[,1],expr1)
rm(expr1)

## Find duplicated genes and calculate the average of them  
dup_ind = which(duplicated(expr[,1]))
dup_ENS = unique(expr[dup_ind,1])

ind = c()
for (i in dup_ENS){
  ind_i = which(expr[,1] == i)
  ave_dup_ENS_i = apply(expr[ind_i,-1],2,mean)
  expr[ind_i[1],-1] = ave_dup_ENS_i
  ind = c(ind,ind_i[-1])
}
expr = expr[-ind,]

#'@Build_expression_matrix_[sample*genes].......................................

rownames(expr) = expr[,1]
expr = expr[,-1]
expr = t(expr)

## gene_ids common between gene_transfer & expression matrix

gene_transfer1 = gene_transfer[,c("gene_id","gene_name")]
gene_transfer1 = gene_transfer1[!duplicated(gene_transfer1[,1]),]

intersect_gene_id = intersect(gene_transfer1$gene_id, colnames(expr))
expr = expr[,intersect_gene_id]
colnames(expr) = gene_transfer1[gene_transfer1$gene_id %in% intersect_gene_id,2]


# Expression matrix normalization
### z_score normalization
expr_norm = scale(expr)
### min-max normalization
#expr_norm_2 = apply(expr,2, function(x){return((x-min(x))/(max(x)-min(x)))})

#'@Build_response_matrix_[sample*Drug]..........................................

cell_id = rownames(expr_norm)
drug_name = unique(response$name)
#sen = matrix(0,length(cell_id),length(drug_name))
#rownames(sen) = cell_id
#colnames(sen) = drug_name
#response = response[!is.na(response$depmap_id),]

### Because for some of the cell lines we do not have the response against 
          #some drugs we put 'NA' for them in line 86-87
# c = 0
# for (i in cell_id){
#   c = c+1
#   print(c)
#   for (j in drug_name){
#     cell_i_drug_j = response$depmap_id == i & response$name == j
#     if (sum(cell_i_drug_j) == 0)
#       sen[i,j] = NA
#     
#     else if (sum(cell_i_drug_j) == 1)
#              sen[i,j] = response[cell_i_drug_j,"auc"]
#     
#     else
#       sen[i,j] = mean(response[cell_i_drug_j,"auc"])
#     
#     }
# }
#saveRDS(sen, file = "Processed_Data/sensitivity_matrix.rds")
sen = readRDS("Data/Processed_Data/sensitivity_matrix.rds")
## Sensitivity matrix normalization
#sen_norm = scale(sen)
#mean_sen_norm = apply(sen_norm,2,function(x){mean(x,na.rm=TRUE)})
#hist(mean_sen_norm)

#### Visualization; just to check the distribution of means and standard deviation across samples and drugs
dist_mean_sample = apply(sen,1,function(x){mean(x,na.rm =TRUE)})
hist(dist_mean_sample)

dist_sd_sample = apply(sen,1,function(x){sd(x,na.rm =TRUE)})
hist(dist_sd_sample)

dist_mean_drug = apply(sen,2,function(x){mean(x,na.rm =TRUE)})
hist(dist_mean_drug)

dist_sd_drug = apply(sen,2,function(x){sd(x,na.rm =TRUE)})
hist(dist_sd_drug)

# Save Data ---------------------------------------------------------------
## 1) sen matrix is saved in line 115 without normalization
saveRDS(expr, file = "Data/Processed_Data/expresion_matrix.rds")
write.table(expr, file = "Data/Processed_Data/expresion_matrix.csv",
            row.names = TRUE, col.names = TRUE, quote = FALSE, sep = ",")
saveRDS(expr_norm, file = "Data/Processed_Data/expresion_normalized_matrix.rds")
saveRDS(gene_transfer1, file = "Data/Processed_Data/gene_transfer.rds")
w = read.csv("Data/Processed_Data/expresion_matrix.csv")
