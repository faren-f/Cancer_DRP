rm(list=ls())
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

# Library -----------------------------------------------------------------
library(igraph)

# Read_Data ---------------------------------------------------------------
GE = readRDS("Processed_Data/S1/expresion_matrix.rds")
#GE_norm = readRDS("Data/Processed_Data/expresion_normalized_matrix.rds")

drug_sensitivity = readRDS("Processed_Data/S1/sensitivity_matrix.rds")

# write.table(drug_sensitivity, file = "Data/Processed_Data_For_Python/Create_cellline_drug_net/Regression/sensitivity_matrix.csv",
#             row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")

#sen = read.csv("sensitivity_matrix.csv")
# Similarity calculation ------------------------------------------------------------

#'@Finding_High_Varuated_genes

# var_GE = apply(GE,2,var)
# var_GE = sort(var_GE,decreasing = TRUE)
# hist(var_GE,100)
# thr_var_GE = 0.8
# abline(v = thr_var_GE, col = "red")
# GE = GE[,var_GE > thr_var_GE]
# write.table(GE, file = "Data/Processed_Data_For_Python/Create_cellline_drug_net/Regression/gene_expresion_highVar.csv",
#             row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")

# Normalization Gene Expresion
#GE = apply(GE,2, function(x){return((x-min(x))/(max(x)-min(x)))})
#GE = scale(GE)

#'@Cell-line_Cell-line_Similarity

sim_GE = abs(cor(t(GE)))

## Distribution of cell line-cell line similarity matrix
hist(sim_GE,30)
thr_binarize = 0.15
abline(v =thr_binarize, col = "red")

## Binarizing the cell line-cell line similarity matrix
sim_GE[sim_GE>thr_binarize] = 1
sim_GE[sim_GE<thr_binarize] = 0

## Assign Node numbers to celllines 
cellline_Node = cbind(rownames(sim_GE),
                      c(0:(nrow(sim_GE)-1)))

rownames(sim_GE) = cellline_Node[,2]
colnames(sim_GE) = cellline_Node[,2]

## Make graph from adj to build edge_index
Graph_sample_sim = graph_from_adjacency_matrix(sim_GE,
                                               mode = "undirected",diag = FALSE)
edge_index_cellline_1 = t(as_edgelist(Graph_sample_sim, names = TRUE))
#'@edge_index_cellline
edge_index_cellline = cbind(edge_index_cellline_1,edge_index_cellline_1[c(2,1),])

#'@node_attr_cellline
node_attr_cellline = GE

#'@Finding_binarized_FingerPrints(node_attr_drug)

Fingerprints = readRDS("Processed_Data/S2/Fingerprints.rds")  

FP = matrix(0, length(Fingerprints), 1024)
for (i in 1:length(Fingerprints)){
  FP_bits_on = Fingerprints[[i]]
  FP[i,FP_bits_on@bits] = 1
}

node_attr_drug = FP

#'@Drug_Drug_Similarity

## Distribution of drug-drug similarity matrix
sim_FP = readRDS("Processed_Data/S2/Fingerprints_sim.rds") 

hist(sim_FP,150)
thr_binarize = 0.12
abline(v = thr_binarize, col = "red")

## Binarizing the drug-drug similarity matrix
sim_FP[sim_FP>thr_binarize] = 1
sim_FP[sim_FP<thr_binarize] = 0

## Assign Node numbers to drugs 
drug_Node = cbind(rownames(sim_FP), c(0:(nrow(sim_FP)-1)))

rownames(sim_FP) = drug_Node[,2]
colnames(sim_FP) = drug_Node[,2]

## Make graph from adj to build edge_index_drug
Graph_drug_sim = graph_from_adjacency_matrix(sim_FP,
                                             mode = "undirected",diag = FALSE)
edge_index_drug_1 = t(as_edgelist(Graph_drug_sim, names = TRUE))

#'@edge_index_drug
edge_index_drug = cbind(edge_index_drug_1,edge_index_drug_1[c(2,1),])

#'@node_attr_drug
node_attr_drug = FP


#'@edge_index_cellline_drug

sen = drug_sensitivity

## Assign Node numbers to sen $ sen_na matrices

rownames(sen) = cellline_Node[,2]
colnames(sen) = drug_Node[,2]

#'@edge_index_cellline_drug
edge_index_cellline_drug = c()
edge_lable_cellline_drug = c()

for (i in 1:nrow(sen)){
  ind1 = names(which(!is.na(sen[i,])))
  edge_i_res = rbind(cellline = rep(i-1, length(ind1)), drug = ind1)
  edge_index_cellline_drug = cbind(edge_index_cellline_drug,edge_i_res) 
  edge_lable_cellline_drug = c(edge_lable_cellline_drug,sen[i,ind1])
}

# Normalization
#edge_lable_cellline_drug = scale(edge_lable_cellline_drug)
#edge_lable_cellline_drug = (edge_lable_cellline_drug-min(edge_lable_cellline_drug))/
  #(max(edge_lable_cellline_drug)-min(edge_lable_cellline_drug))
#hist(edge_lable_cellline_drug,30)
# Save_Data ---------------------------------------------------------------

#'@Save_data_cell_line
write.table(node_attr_cellline, file = "Processed_Data_For_Python/Create_cellline_drug_net/Regression/node_attr_cellline.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
write.table(edge_index_cellline, file = "Processed_Data_For_Python/Create_cellline_drug_net/Regression/edge_index_cellline.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")


#'@Save_data_drug
write.table(node_attr_drug, file = "Processed_Data_For_Python/Create_cellline_drug_net/Regression/node_attr_drug.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
write.table(edge_index_drug, file = "Processed_Data_For_Python/Create_cellline_drug_net/Regression/edge_index_drug.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")


#'@Save_data_cellline_drug
write.table(edge_index_cellline_drug, file = "Processed_Data_For_Python/Create_cellline_drug_net/Regression/edge_index_cellline_drug.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
write.table(edge_lable_cellline_drug, file = "Processed_Data_For_Python/Create_cellline_drug_net/Regression/edge_lable_cellline_drug.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")

