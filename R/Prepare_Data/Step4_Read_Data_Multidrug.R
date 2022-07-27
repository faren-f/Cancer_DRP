rm(list=ls())
setwd("~/Desktop/Codes/Cancer_DRP/R/")

# Library -----------------------------------------------------------------
library(igraph)

# Read_Data ---------------------------------------------------------------
GE = readRDS("Data/Processed_Data/expresion_matrix.rds")
#GE_norm = readRDS("Data/Processed_Data/expresion_normalized_matrix.rds")

drug_sensitivity = readRDS("Data/Processed_Data/sensitivity_matrix.rds")


# Similarity calculation ------------------------------------------------------------

#'@Cell-line_Cell-line_Similarity

var_GE = apply(GE,2,var)
var_GE = sort(var_GE,decreasing = TRUE)
hist(var_GE,100)
thr_var_GE = 0.7
abline(v = thr_var_GE, col = "red")

GE = GE[,var_GE > thr_var_GE]
sim_GE = cor(t(GE))

## Distribution of cell line-cell line similarity matrix
hist(sim_GE,30)
thr_binarize = 0.68
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



#'@Drug_Drug_Similarity
Fingerprints = readRDS("Data/Processed_Data/Fingerprints.rds")  

FP = matrix(0, length(Fingerprints), 1024)
for (i in 1:length(Fingerprints)){
  FP_bits_on = Fingerprints[[i]]
  FP[i,FP_bits_on@bits] = 1
}

node_attr_drug = FP
## Distribution of drug-drug similarity matrix
sim_FP = readRDS("Data/Processed_Data/Fingerprints_sim.rds") 

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

## Distribution of drug_sensitivity matrix
hist(drug_sensitivity,100)
thr_sen = 1
abline(v = thr_sen, col = "red")

## Binarizing the drug_sensitivity matrix
sen = drug_sensitivity
sen[sen>1] = 1
sen[sen<1] = 0
sen[is.na(sen)] = 0

sen_na = drug_sensitivity
sen_na[sen_na>1] = 1
sen_na[sen_na<1] = -1
sen_na[is.na(sen_na)] = 0

## Assign Node numbers to sen $ sen_na matrices

rownames(sen) = cellline_Node[,2]
rownames(sen_na) = cellline_Node[,2]

colnames(sen) = drug_Node[,2]
colnames(sen_na) = drug_Node[,2]

#'@edge_index_cellline_drug
edge_index_cellline_drug = c()
edge_index_cellline_drug_NA = c()

for (i in 1:nrow(sen)){
  
  ind1 = names(which(sen[i,] == 1))
  edge_i_res = rbind(cellline = rep(i-1, length(ind1)), drug = ind1)
  edge_index_cellline_drug = cbind(edge_index_cellline_drug,edge_i_res)  
  
  ind2 = names(which(sen_na[i,] == 0))
  edge_i_na = rbind(cellline = rep(i-1, length(ind2)), drug = ind2)
  edge_index_cellline_drug_NA = cbind(edge_index_cellline_drug_NA,edge_i_na)  
}


# Save_Data ---------------------------------------------------------------

#'@Save_data_cell_line
write.table(node_attr_cellline, file = "Data/Processed_Data_For_Python/Create_cellline_drug_net/Classifier/node_attr_cellline.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
write.table(edge_index_cellline, file = "Data/Processed_Data_For_Python/Create_cellline_drug_net/Classifier/edge_index_cellline.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")


#'@Save_data_drug
write.table(node_attr_drug, file = "Data/Processed_Data_For_Python/Create_cellline_drug_net/Classifier/node_attr_drug.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
write.table(edge_index_drug, file = "Data/Processed_Data_For_Python/Create_cellline_drug_net/Classifier/edge_index_drug.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")


#'@Save_data_cellline_drug
write.table(edge_index_cellline_drug, file = "Data/Processed_Data_For_Python/Create_cellline_drug_net/Classifier/edge_index_cellline_drug.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
write.table(edge_index_cellline_drug_NA, file = "Data/Processed_Data_For_Python/Create_cellline_drug_net/Classifier/edge_index_cellline_drug_NA.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")

