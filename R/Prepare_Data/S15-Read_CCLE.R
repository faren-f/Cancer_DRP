#                   Created on Thu Jan 02 10:16 2021

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This script reads raw RNAseq data and drug sensitivity data from CCLE and convert
# gene-ids to gene symbols and do log2 normalization on RNAseq data. finally it prepare drug 
# sensitivity matrix and gene expression matrix.

rm(list=ls())
library('rtracklayer')
library(dplyr)
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
response = read.csv("Raw_data/CCLE/Drug_Sensitivity/CCLE_NP24.2009_Drug_data_2015.02.24.csv")

RNAseq = read.table("Raw_data/CCLE/Expression/RNA_seq/CCLE_RNAseq_rsem_transcripts_tpm_20180929.txt.gz",
                     header = TRUE, check.names = FALSE)
gene_transfer = import("Raw_data/CCLE/Expression/RNA_seq/gencode.v19.genes.v7_model.patched_contigs.gtf.gz")
gene_transfer = data.frame(gene_transfer)
cellline_info = read.csv("Raw_data/CCLE/Expression/RNA_seq/sample_info.csv")


# Pre-processing ----------------------------------------------------------
## Log normalization of genes
expr_raw = RNAseq[,c(-1,-2)]
expr = log2(expr_raw + 1)


## Remove cell lines that do not exist in response from expression
ccle_name_intersect = colnames(expr) %in% response$CCLE.Cell.Line.Name
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
depid_name = cellline_info[,c(1,4)]
depid_name = mutate_all(depid_name, na_if,"")
depid_name = depid_name[!duplicated(depid_name[,2]),]
depid_name = depid_name[!is.na(depid_name[,2]),]
rownames(depid_name) = depid_name[,2]

i = intersect(colnames(expr1),depid_name[,2])

depid_name = depid_name[i,]
expr1 = expr1[,i]
colnames(expr1) = depid_name[,1] 

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

intersect_gene_id = intersect(colnames(Expr), gene_transfer1$gene_id)
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


ccle_name_Expr = depid_name[,2]
drug_name = unique(response$Compound)
IC50 = matrix(0,length(ccle_name_Expr),length(drug_name))
AA = matrix(0,length(ccle_name_Expr),length(drug_name))

rownames(IC50) = ccle_name_Expr
colnames(IC50) = drug_name

rownames(AA) = ccle_name_Expr
colnames(AA) = drug_name

### Because for some of the cell lines we do not have the response against 
# some drugs we put 'NA' for them in line 86-87
c = 0
for (i in ccle_name_Expr){
  c = c+1
  print(c)
  
  for (j in drug_name){
    cell_i_drug_j = response$CCLE.Cell.Line.Name == i & response$Compound == j
    
    if (sum(cell_i_drug_j) == 0)
      IC50[i,j] = NA
      AA[i,j] = NA

    else if (sum(cell_i_drug_j) == 1){
             IC50[i,j] = response[cell_i_drug_j,11]
             AA[i,j] = response[cell_i_drug_j,13]
    
    }else{
      IC50[i,j] = mean(response[cell_i_drug_j,11])
      AA[i,j] = mean(response[cell_i_drug_j,13])
      }
   }
}
rownames(IC50) = rownames(Expr)
rownames(AA) = rownames(Expr)


saveRDS(IC50, file = "Processed_Data/S15/sensitivity_matrix_IC50.rds")
saveRDS(AA, file = "Processed_Data/S15/sensitivity_matrix_Activity_Area.rds")

#IC50 = readRDS("Processed_Data/S15/sensitivity_matrix_IC50.rds")
#AA = readRDS("Processed_Data/S15/sensitivity_matrix_Activity_Area.rds")

#### Visualization; just to check the distribution of means and standard deviation across samples and drugs
# dist_mean_sample = apply(IC50,1,function(x){mean(x,na.rm =TRUE)})
# hist(dist_mean_sample)
#
# dist_sd_sample = apply(IC50,1,function(x){sd(x,na.rm =TRUE)})
# hist(dist_sd_sample)
#
# dist_mean_drug = apply(IC50,2,function(x){mean(x,na.rm =TRUE)})
# hist(dist_mean_drug)
#
# dist_sd_drug = apply(IC50,2,function(x){sd(x,na.rm =TRUE)})
# hist(dist_sd_drug)

# Save Data ---------------------------------------------------------------
## 1) IC50, activity area matrices are saved in line 178-179
saveRDS(Expr, file = "Processed_Data/S15/expresion_matrix.rds")
write.table(Expr, file = "Processed_Data/S15/expresion_matrix.csv",
            row.names = TRUE, col.names = TRUE, quote = FALSE, sep = ",")

saveRDS(gene_transfer1, file = "Processed_Data/S15/gene_transfer.rds")

