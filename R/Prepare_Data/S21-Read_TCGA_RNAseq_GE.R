rm(list=ls())
library("readxl")

# Read clinical data from table (2016-Ding)
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
N_Patients = readRDS("Processed_data/S21/Num_Patients_in_each_cancer.rds")
Response_table = read_excel("Raw_data/TCGA/2016-Ding/bioinfo16_supplementary_tables.xlsx",
                      sheet = 3,na = "---")
colnames(Response_table) = Response_table[2,]
Response_table = Response_table[c(-1,-2,-3),]
Response_table$drug_name = tolower(Response_table$drug_name) 
Response_table = data.frame(Response_table)

Cancer_types = unique(Response_table$Cancer)
Cancer = sapply(strsplit(Cancer_types,"\\("), FUN = function(x){return(x[2])})
Cancer = sapply(strsplit(Cancer,"\\)"), FUN = function(x){return(x[1])})
S = c()
for(i in Cancer_types){
  S = c(S,sum(Response_table$Cancer == i))
}
x = data.frame(c = Cancer,s = S)
Cancer_type = data.frame(lapply(x, rep, x$s))
Response_table$Cancer_type = Cancer_type[,1]

## Read RNAseq data from (http://firebrowse.org) 
BRCA = read.table("Raw_data/TCGA/RNAseq/BRCA/BRCA.rnaseqv2_RSEM_genes_normalized.txt",
                      fill = TRUE, header=FALSE)
BLCA = read.table("Raw_data/TCGA/RNAseq/BLCA/BLCA.rnaseqv2_RSEM_genes_normalized.txt",
                      fill = TRUE, header=FALSE)

Cancers = c("BRCA","BLCA")
TCGA_GE = c()
TT = c()
for (i in Cancers){
  Cancer_i = get(i)
  Cancer_i = Cancer_i[-2,]
  patient_all = Cancer_i[1,]
  patient_all = substr(Cancer_i[1,],1,12)
  Cancer_i[1,] = patient_all
  
  patient_with_response = Response_table[which(Response_table$Cancer_type == i),2]
  Intersect = intersect(patient_all, patient_with_response)
  
  patient_with_response = Response_table[Response_table$bcr_patient_barcode %in% Intersect,]
  Cancer_i = Cancer_i[,c(1,which(patient_all %in% Intersect))]
  dup = which(duplicated(as.character(Cancer_i[1,])))
  dup = dup-1
  
  colnames(Cancer_i) = Cancer_i[1,]
  Cancer_i = Cancer_i[-1,]
  rownames(Cancer_i) = Cancer_i[,1]
  Cancer_i = Cancer_i[,c(-1)]
  Cancer_i = t(Cancer_i)
  #B = data.frame((rownames(Cancer_i)))
  Cancer_i_2 = apply(Cancer_i,2,as.numeric)
  
  for (k in dup){
    rep = apply(Cancer_i_2[c(k-1,k),],2,mean)
    Cancer_i_2[k-1,] = rep
  }
  
  Cancer_i_2 = Cancer_i_2[-dup,]
  Cancer_i = Cancer_i[-dup,]
  rownames(Cancer_i_2) = rownames(Cancer_i)
  colnames(Cancer_i_2) = colnames(Cancer_i)
  Cancer_i = Cancer_i_2
  
  col = colnames(Cancer_i)
  col = sapply(strsplit(col,"\\|"),FUN = function(x){return(x[1])})
  colnames(Cancer_i) = col
  Cancer_i = log2(Cancer_i + 1)
  
  GE = readRDS("Processed_data/S1/expresion_matrix.rds")
  intersect_genes = intersect(colnames(GE),colnames(Cancer_i))
  
  GE = GE[,intersect_genes]
  Cancer_i = Cancer_i[,intersect_genes]
  TT = c(TT,rep(i,nrow(Cancer_i)))
  TCGA_GE = rbind(TCGA_GE,Cancer_i)
}
  
#hist(Cancer_i,100)



