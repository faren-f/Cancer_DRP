rm(list=ls())
library(PharmacoGx)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
sen = read.table("Raw_data/GDSC1/Drug_Response/PANCANCER_IC_Mon Nov 14 10_28_56 2022.csv",
                  fill = TRUE, sep = ",",header = TRUE)

GDSC = readRDS("Processed_data/S38/GDSC_PharmacoGx.rds")

CellName = phenoInfo(GDSC,"rna") 
GeneName = featureInfo(GDSC, "rna")                               
GE = t(molecularProfiles(GDSC, "rna"))                         
rownames(GE) = CellName$Characteristics.cell.line.
colnames(GE) = GeneName$Symbol
#hist(GE)

#DrugInfo = drugInfo(GDSC)
#CellInfo = cellInfo(GDSC)                                         
#CellInfo = CellInfo[CellName$Characteristics.cell.line.,]        


## Remove cell lines that do not exist in response from expression
gdsc_name_intersect = rownames(GE) %in% sen$Cell.Line.Name
GE = GE[gdsc_name_intersect,]
dim(GE)



#'@Build_response_matrix_[sample*Drug]..........................................

cellline_name = rownames(GE)
drug_name = unique(sen$Drug.Name)
AUC = matrix(0,length(cellline_name),length(drug_name))
IC50 = matrix(0,length(cellline_name),length(drug_name))

rownames(AUC) = cellline_name
colnames(AUC) = drug_name

rownames(IC50) = cellline_name
colnames(IC50) = drug_name

## Because for some of the cell lines we do not have the response against
# some drugs we put 'NA' for them
c = 0
for (i in cellline_name){
  c = c+1
  print(c)
  
  for (j in drug_name){
    cell_i_drug_j = sen$Cell.Line.Name == i & sen$Drug.Name == j
    
    if (sum(cell_i_drug_j) == 0){
      AUC[i,j] = NA
      IC50[i,j] = NA
      
    }else if (sum(cell_i_drug_j) == 1){
      AUC[i,j] = sen[cell_i_drug_j,"AUC"]
      IC50[i,j] = sen[cell_i_drug_j,"IC50"]

    }else{
      AUC[i,j] = mean(sen[cell_i_drug_j,"AUC"])
      IC50[i,j] = mean(sen[cell_i_drug_j,"IC50"])

    }
  }
}
#saveRDS(AUC, file = "Processed_Data/S41/GDSC_sensitivity_matrix_AUC.rds")
#saveRDS(IC50, file = "Processed_Data/S41/GDSC_sensitivity_matrix_IC50.rds")
#saveRDS(GE, file = "Processed_Data/S41/GDSC_expresion_matrix.rds")

