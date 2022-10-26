rm(list=ls())
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

# Library -----------------------------------------------------------------
library(igraph)

# Read_Data ---------------------------------------------------------------
#GE = readRDS("Processed_Data/S1/expresion_matrix.rds")
GE = readRDS("Processed_data/Other/GE_PRISM_CommonGeneswith_TCGA.rds")
l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
GE = GE[,colnames(GE)%in%l1000_genes]
drug_sensitivity = readRDS("Processed_Data/S1/sensitivity_matrix.rds")
#drug_sensitivity = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")

#write.table(drug_sensitivity, file = "Data/Processed_Data_For_Python/Create_cellline_drug_net/Regression/sensitivity_matrix.csv",
#           row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")

# Similarity calculation ------------------------------------------------------------

#'@Normalization_of_GeneExpresion

#GE = apply(GE,2, function(x){return((x-min(x))/(max(x)-min(x)))})
GE = scale(GE)

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
cellline_Node = cbind(rownames(sim_GE),c(0:(nrow(sim_GE)-1)))

rownames(sim_GE) = cellline_Node[,2]
colnames(sim_GE) = cellline_Node[,2]

## Make graph from adj to build edge_index
Graph_sample_sim = graph_from_adjacency_matrix(sim_GE,mode = "undirected",diag = FALSE)

edge_index_cellline_1 = t(as_edgelist(Graph_sample_sim, names = TRUE))
#'@edge_index_cellline
edge_index_cellline = cbind(edge_index_cellline_1,edge_index_cellline_1[c(2,1),])

#'@node_attr_cellline
node_attr_cellline = GE

#'@Finding_binarized_FingerPrints(node_attr_drug)

Fingerprints = readRDS("Processed_Data/S2/Fingerprints.rds")  
good_drugs = readRDS("Processed_data/Other//good_drugs_in_PRISM.rds")
Fingerprints = Fingerprints[good_drugs[,1]]

FP = matrix(0, length(Fingerprints), 1024)
for (i in 1:length(Fingerprints)){
  FP_bits_on = Fingerprints[[i]]
  FP[i,FP_bits_on@bits] = 1
}

node_attr_drug = FP

#'@Drug_Drug_Similarity

drug_sensitivity = drug_sensitivity[,rownames(good_drugs)]

Fingerprints_sim = fingerprint::fp.sim.matrix(Fingerprints, method='tanimoto')
rownames (Fingerprints_sim) =  rownames(good_drugs)
colnames (Fingerprints_sim) =  rownames(good_drugs)

#saveRDS(Fingerprints_sim,"Processed_Data/S37/Fingerprints_sim.rds")  


## Distribution of drug-drug similarity matrix
sim_FP = Fingerprints_sim

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

sensitivity = apply(drug_sensitivity,2, function(x){
  m = mean(x, na.rm = TRUE)
  s = sd(x, na.rm = TRUE)
  x = (x-m)/s
  return(x)})

#sensitivity = t(sensitivity)
sen = sensitivity

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
hist(edge_lable_cellline_drug)
edge_lable = edge_lable_cellline_drug
q1 = quantile(edge_lable,prob=0.05)
q3 = quantile(edge_lable,prob=0.95)

edge_lable = (edge_lable-q1)/(q3-q1)
#edge_lable_cellline_drug = edge_lable
#hist(edge_lable, 100)
edge_lable[edge_lable>1] = 1
edge_lable[edge_lable<0] = 0

#hist(edge_lable, 100)

edge_weight_cellline_drug = 1-edge_lable
#edge_weight_cellline_drug = rep(1, length(edge_lable))
#hist(edge_weight_cellline_drug, 100)

# Normalization
#edge_lable_cellline_drug = scale(edge_lable_cellline_drug)
#edge_lable_cellline_drug = (edge_lable_cellline_drug-min(edge_lable_cellline_drug))/
#(max(edge_lable_cellline_drug)-min(edge_lable_cellline_drug))
#hist(edge_lable_cellline_drug,30)

# Celllines-Drugs that are NA
Na_sen = ifelse(apply(sensitivity,2,function(x){return(is.na(x))}),0,1)

# Save_Data ---------------------------------------------------------------

#'@Save_data_cell_line
write.table(node_attr_cellline, file = "Processed_Data/S37-P/node_attr_cellline.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
write.table(edge_index_cellline, file = "Processed_Data/S37-P/edge_index_cellline.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")


#'@Save_data_drug
write.table(node_attr_drug, file = "Processed_Data/S37-P/node_attr_drug.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
write.table(edge_index_drug, file = "Processed_Data/S37-P/edge_index_drug.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")


#'@Save_data_cellline_drug
write.table(edge_index_cellline_drug, file = "Processed_Data/S37-P/edge_index_cellline_drug.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
write.table(edge_lable_cellline_drug, file = "Processed_Data/S37-P/edge_lable_cellline_drug.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
write.table(edge_weight_cellline_drug, file = "Processed_Data/S37-P/edge_weight_cellline_drug.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")

write.table(Na_sen, file = "Processed_Data/S37-P/Na_sen.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")

write.table(sensitivity, file = "Processed_Data/S37-P/drug_sensitivity.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")




