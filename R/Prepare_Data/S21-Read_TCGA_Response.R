rm(list=ls())
library("readxl")

# Read TCGA clinical data from table (2016-Ding)
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
Response_table = read_excel("Raw_data/TCGA/2016-Ding/bioinfo16_supplementary_tables.xlsx",
                      sheet = 3,na = "---")
colnames(Response_table) = Response_table[2,]
Response_table = Response_table[c(-1,-2,-3),]
Response_table$drug_name = tolower(Response_table$drug_name) 
Response_table = data.frame(Response_table)

# Drug names that are available in TCGA
drugs = unique(Response_table$drug_name)

# drugs that are common between TCGA and PRISM
sen = readRDS("Processed_data/S1/sensitivity_matrix_AUC.rds")
TCGA_PRISM_drugs = intersect(drugs,colnames(sen))
saveRDS(TCGA_PRISM_drugs,"Processed_data/S21/Drugs_TCGA@PRISM.rds")

# Patients & Cancer type in TCGA
Cancer_Patient = Response_table[!duplicated(Response_table$bcr_patient_barcode),1:2]

# Drug response measures
res_measure = unique(Response_table$measure_of_response)
res_binarized = rep(0, nrow(Response_table))

for (i in 1:length(res_measure)){
    res_binarized[which(Response_table$measure_of_response == res_measure[i])] = i
}
Response_table$response_binarized = res_binarized

# Patients-Response matrix
response_mat = matrix(0,nrow(Cancer_Patient),length(TCGA_PRISM_drugs))
rownames(response_mat) = Cancer_Patient[,2]
colnames(response_mat) = TCGA_PRISM_drugs

for(i in Cancer_Patient[,2]){
  for (j in TCGA_PRISM_drugs){
    I = Response_table$bcr_patient_barcode == i & Response_table$drug_name == j
    if(sum(I)==1){
      response_mat[i,j] = Response_table[I,15]
      
    }else if(sum(I)==2){
      if(Response_table[which(I)[1],15]==Response_table[which(I)[2],15])
        response_mat[i,j] = Response_table[which(I)[1],15]
    }else{
      response_mat[i,j] = NA
    }
  }
}
saveRDS(response_mat,"Processed_data/S21/Drug_Response_matrix_TCGA.rds")

# Cancers types in TCGA that have drug responses
cancers = unique(Response_table$Cancer)
saveRDS(cancers,"Processed_data/S21/Cancer_types_TCGA.rds")
# Number of patients in each cancer
N_P = c()
for (i in cancers){
  N_P_i = Response_table[Response_table$Cancer==i,2]
  N_P_i = N_P_i[!is.na(N_P_i)]
  N_P = c(N_P,length(unique(N_P_i)))
}
saveRDS(N_P,"Processed_data/S21/Num_Patients_in_each_cancer.rds")

