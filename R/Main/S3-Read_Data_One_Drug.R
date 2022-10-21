rm(list=ls())
setwd("~/Desktop/Codes/Cancer_DRP/R/")


# Library -----------------------------------------------------------------
library(igraph)

# Read_Data ---------------------------------------------------------------
expr = readRDS("Data/Processed_Data/expresion_matrix.rds")
expr_norm = readRDS("Data/Processed_Data/expresion_normalized_matrix.rds")

sen = readRDS("Data/Processed_Data/sensitivity_matrix.rds")
#sample_sim_exp = readRDS("Data/Processed_Data/Sample_Sim_Exp.rds")


## Find drug responses for drug i and remove NA samples 

#max_cor = c()                             ##from sen and expr matrix 
#for(i in 1:ncol(sen)){
  #sen_drug_i_with_NA = sen[,i]
  #indeces_without_NA = which(!is.na(sen_drug_i_with_NA))
  #sen_drug_i = sen_drug_i_with_NA[indeces_without_NA]
  #expr_i = expr[indeces_without_NA,]
  
  #cor_Xy = apply(expr_i,2,function(x){abs(cor(x,as.numeric(sen_drug_i)))})
  #max_cor = c(max_cor, max(cor_Xy))
#}
#max_cor = data.frame(max_cor)
#saveRDS(max_cor,"Data/Processed_Data/max_cor_genes_auc.rds")
max_cor = readRDS("Data/Processed_Data/max_cor_genes_auc.rds")
sen_drug_i = sen[,60]
indeces_without_NA = which(!is.na(sen_drug_i))
sen_drug_i = sen_drug_i[indeces_without_NA]
expr_i = expr[indeces_without_NA,]

#'@Output
y = sen_drug_i
hist(y,30)
thr_y = 0 
abline(v = thr_y,col = "red")
#y = ifelse(y > thr_y,1,0) # For Classification

expr_i = expr[indeces_without_NA,]

### Feature selection
#with corr
thr_cor = 0.15
cor_Xy = apply(expr_i,2,function(x){abs(cor(x,as.numeric(y)))})
hist(cor_Xy,20)
ind_high_cor = which(cor_Xy>thr_cor)
expr_i = expr_i[,ind_high_cor]

#with sd
#sd_expr = apply(expr_i, 2, sd)
#thr_feat = 1.2
#hist(sd_expr,30)
#abline(v = thr_feat,col = "red")
#ind_high_sd = which(sd_expr>thr_feat)

#expr_i = expr_i[,ind_high_sd]

## Similarity calculation
sim_expr_i = cor(t(expr_i))


## Distribution of cell line-cell line similarity matrix
hist(sim_expr_i,30)
thr_binarize = 0.6
abline(v =thr_binarize, col = "red")

## Binarizing the cell line-cell line similarity matrix
sim_expr_i[sim_expr_i>thr_binarize] = 1
sim_expr_i[sim_expr_i<thr_binarize] = 0

## Assign Node numbers to samples 
Sample_Node = cbind(rownames(sim_expr_i),
                    c(0:(nrow(sim_expr_i)-1)))

rownames(sim_expr_i) = Sample_Node[,2]
colnames(sim_expr_i) = Sample_Node[,2]

## Make graph from adj to build edge_index
Graph_sample_sim = graph_from_adjacency_matrix(sim_expr_i,
                                               mode = "undirected",diag = FALSE)
edge_index_1 = t(as_edgelist(Graph_sample_sim, names = TRUE))
#'@edge_index
edge_index = cbind(edge_index_1,edge_index_1[c(2,1),])

#'@node_attr
node_attr = expr_i


# Save_Data ---------------------------------------------------------------

write.table(edge_index, file = "Data/Processed_Data_For_Python/Create_cellline_net/edge_index.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")

write.table(node_attr, file = "Data/Processed_Data_For_Python/Create_cellline_net/node_attr.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")

write.table(y, file = "Data/Processed_Data_For_Python/Create_cellline_net/y.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")



