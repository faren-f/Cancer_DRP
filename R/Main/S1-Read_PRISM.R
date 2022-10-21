
#                  Created on Thu Jan 02 10:16 2021

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This script reads raw data from CCLE cell lines and PRISM drug sensitivity and convert
# gene-ids to gene symbols and do log2 normalization on RNAseq data. finally it prepare drug 
# sensitivity matrix and gene expression matrix.

rm(list = ls())
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

# Library -----------------------------------------------------------------
library('rtracklayer')
library(ggplot2)
# Read Data ---------------------------------------------------------------
cellline_info = read.csv("Raw_data/PRISM/Secondary/secondary-screen-cell-line-info.csv")

response = read.csv("Raw_data/PRISM/Secondary/secondary-screen-dose-response-curve-parameters.csv")

RNAseq = read.table("Raw_data/PRISM/Expression/RNA_seq/CCLE_RNAseq_rsem_transcripts_tpm_20180929.txt.gz",
                    header = TRUE, check.names = FALSE)

# ccle = read.csv("Raw_data/PRISM/Expression/RNA_seq/CCLE_expression.csv",
#                 header = FALSE)

#RNAseq = read.table("Raw_data/PRISM/Expression/RNA_seq/CCLE_RNAseq_genes_rpkm_20180929.gct.txt", 
#skip = 2, header = TRUE, sep = "\t")

gene_transfer = import("Raw_data/PRISM/Expression/RNA_seq/gencode.v19.genes.v7_model.patched_contigs.gtf.gz")
gene_transfer = data.frame(gene_transfer)

# Finding drug targets for all drugs
drug_targets = response[!duplicated(response$name),c(12,14)]
rownames(drug_targets) = 1:nrow(drug_targets)
# Pre-processing ----------------------------------------------------------
## Log normalization of genes
expr_raw = RNAseq[,c(-1,-2)]
expr = log2(expr_raw + 1)


## Remove cell lines that do not exist in response from expression
ccle_name_intersect = colnames(expr) %in% response$ccle_name
expr = expr[,ccle_name_intersect]


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
hist(mean_expr,100,xlim = c(0,5))
abline(v = 0.2, col ="red")
expr = expr[mean_expr > 0.2,]

### Depict Figure 2.A
#mean_expr = mean_expr[mean_expr>1]
#pdf(file = "Figures/Fig2_A_expr_mean_filter.pdf", width = 5, height = 3.5)
#hist(mean_expr, 100, main = "Histogram of expression mean",xlab = " ")
#dev.off()

## Remove expressions with low std 
# sd_expr = apply(expr[,-1], 1, sd)
# hist(sd_expr,100)
# abline(v = 0.45, col ="red")
# expr = expr[sd_expr > 0.45,]

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

## Convert ccle_name to depmap_id in expression matrix
depid_name = response[,c(2,3)]
depid_name = depid_name[!duplicated(depid_name[,2]),]
depid_name = depid_name[!is.na(depid_name[,2]),]
rownames(depid_name) = depid_name[,2]
dep_id = depid_name[colnames(expr1),1] 
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
Expr = expr[-ind,]

#'@Build_expression_matrix_[sample*genes].......................................

rownames(Expr) = Expr[,1]
Expr = Expr[,-1]
Expr = t(Expr)

## gene_ids common between gene_transfer & expression matrix

gene_transfer1 = gene_transfer[,c("gene_id","gene_name")]
gene_transfer1 = gene_transfer1[!duplicated(gene_transfer1[,1]),]

intersect_gene_id = intersect(gene_transfer1$gene_id, colnames(Expr))
Expr = Expr[,intersect_gene_id]
colnames(Expr) = gene_transfer1[gene_transfer1$gene_id %in% intersect_gene_id,2]
# Removing repeatative gene-symboles 
col_Exp = colnames(Expr)
dup_gene_symbs = unique(col_Exp[duplicated(colnames(Expr))])

ind_extra = c()
for (k in dup_gene_symbs){
  ind_rep = which(colnames(Expr)==k)
  Expr[,ind_rep] = apply(Expr[,ind_rep],1,mean)
  ind_extra = c(ind_extra,ind_rep[-1])
}
Expr = Expr[,-ind_extra]
#'@Build_response_matrix_[sample*Drug]..........................................

cell_id = rownames(Expr)
drug_name = unique(response$name)
AUC = matrix(0,length(cell_id),length(drug_name))
IC50 = matrix(0,length(cell_id),length(drug_name))
EC50 = matrix(0,length(cell_id),length(drug_name))

rownames(AUC) = cell_id
colnames(AUC) = drug_name

rownames(IC50) = cell_id
colnames(IC50) = drug_name

rownames(EC50) = cell_id
colnames(EC50) = drug_name

response = response[!is.na(response$depmap_id),]

## Because for some of the cell lines we do not have the response against
# some drugs we put 'NA' for them in line 86-87
c = 0
for (i in cell_id){
  c = c+1
  print(c)
  
  for (j in drug_name){
    cell_i_drug_j = response$depmap_id == i & response$name == j
    
    if (sum(cell_i_drug_j) == 0){
      AUC[i,j] = NA
      IC50[i,j] = NA
      EC50[i,j] = NA

    }else if (sum(cell_i_drug_j) == 1){
             AUC[i,j] = response[cell_i_drug_j,"auc"]
             IC50[i,j] = response[cell_i_drug_j,"ic50"]
             EC50[i,j] = response[cell_i_drug_j,"ec50"]
             
    }else{
      AUC[i,j] = mean(response[cell_i_drug_j,"auc"])
      IC50[i,j] = mean(response[cell_i_drug_j,"ic50"])
      EC50[i,j] = mean(response[cell_i_drug_j,"ec50"])
      
      }
   }
}

saveRDS(AUC, file = "Processed_Data/S1/sensitivity_matrix_AUC.rds")
saveRDS(IC50, file = "Processed_Data/S1/sensitivity_matrix_IC50.rds")
saveRDS(EC50, file = "Processed_Data/S1/sensitivity_matrix_EC50.rds")

write.table(AUC, file = "Processed_Data/S1/sensitivity_matrix_AUC.csv",
            row.names = TRUE, col.names = TRUE, quote = FALSE, sep = ",")
write.table(IC50, file = "Processed_Data/S1/sensitivity_matrix_IC50.csv",
            row.names = TRUE, col.names = TRUE, quote = FALSE, sep = ",")
write.table(EC50, file = "Processed_Data/S1/sensitivity_matrix_EC50.csv",
            row.names = TRUE, col.names = TRUE, quote = FALSE, sep = ",")


#AUC = readRDS("Processed_Data/S1/sensitivity_matrix_AUC.rds")
#IC50 = readRDS("Processed_Data/S1/sensitivity_matrix_IC50.rds")
#EC50 = readRDS("Processed_Data/S1/sensitivity_matrix_EC50.rds")

#### Visualization; just to check the distribution of means and standard deviation across samples and drugs
# dist_mean_sample = apply(AUC,1,function(x){mean(x,na.rm =TRUE)})
# hist(dist_mean_sample)
# 
# dist_sd_sample = apply(AUC,1,function(x){sd(x,na.rm =TRUE)})
# hist(dist_sd_sample)
# 
# dist_mean_drug = apply(AUC,2,function(x){mean(x,na.rm =TRUE)})
# hist(dist_mean_drug)
# 
# dist_sd_drug = apply(AUC,2,function(x){sd(x,na.rm =TRUE)})
# hist(dist_sd_drug)

# Save Data ---------------------------------------------------------------
## 1) AUC and IC50 matrices are saved in line 179-180
saveRDS(Expr, file = "Processed_Data/S1/expresion_matrix.rds")
write.table(Expr, file = "Processed_Data/S1/expresion_matrix.csv",
            row.names = TRUE, col.names = TRUE, quote = FALSE, sep = ",")

saveRDS(gene_transfer1, file = "Processed_Data/S1/gene_transfer.rds")
saveRDS(drug_targets, file = "Processed_Data/S1/drug_targets.rds")






