rm(list=ls())
setwd("~/Desktop/Codes/Cancer_DRP/R/")

# Library -----------------------------------------------------------------
library(igraph)

# Read_Data ---------------------------------------------------------------
GE = readRDS("Data/Processed_Data/expresion_matrix.rds")
#GE_norm = readRDS("Data/Processed_Data/expresion_normalized_matrix.rds")
drug_sensitivity = readRDS("Data/Processed_Data/sensitivity_matrix.rds")


# Similarity calculation ------------------------------------------------------------

#'@Finding_High_Varuated_genes

var_GE = apply(GE,2,var)
var_GE = sort(var_GE,decreasing = TRUE)
hist(var_GE,100)
thr_var_GE = 0.8
abline(v = thr_var_GE, col = "red")

GE = GE[,var_GE > thr_var_GE]

# Normalization Gene Expresion
GE = scale(GE)
#GE = apply(GE,2, function(x){return((x-min(x))/(max(x)-min(x)))})


#'@Finding_binarized_FingerPrints
Fingerprints = readRDS("Data/Processed_Data/Fingerprints.rds")  

FP = matrix(0, length(Fingerprints), 1024)
for (i in 1:length(Fingerprints)){
  FP_bits_on = Fingerprints[[i]]
  FP[i,FP_bits_on@bits] = 1
}


#'@cellline_drug

## Assign index to cell line-drug pair
sensitivity = drug_sensitivity
cellline_Node = cbind(rownames(GE), c(0:(nrow(GE)-1)))
drug_Node = cbind(colnames(sensitivity), c(0:(ncol(sensitivity)-1)))

rownames(sensitivity) = cellline_Node[,2]
colnames(sensitivity) = drug_Node[,2]

#'@cellline_drug
cellline_drug_index = c()
sen = c()

for (i in 1:nrow(sensitivity)){
  ind1 = names(which(!is.na(sensitivity[i,])))
  edge_i_res = rbind(cellline = rep(i-1, length(ind1)), drug = ind1)
  cellline_drug_index = cbind(cellline_drug_index,edge_i_res) 
  sen = c(sen,sensitivity[i,ind1])
}

# Normalization
sen = scale(sen)
#sen = apply(sen,2, function(x){return((x-min(x))/(max(x)-min(x)))})
hist(sen,30)


# Save_Data ---------------------------------------------------------------

#'@Save_data_cell_line
write.table(GE, file = "Data/Processed_Data_For_Python/Create_Data/ENet/GE.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
write.table(FP, file = "Data/Processed_Data_For_Python/Create_Data/ENet/FP.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
write.table(sen, file = "Data/Processed_Data_For_Python/Create_Data/ENet/sen.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
write.table(cellline_drug_index, file = "Data/Processed_Data_For_Python/Create_Data/ENet/cellline_drug_index.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")

