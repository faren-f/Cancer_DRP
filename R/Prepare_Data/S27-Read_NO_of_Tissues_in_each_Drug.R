rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
TCGA_Patients = readRDS("Processed_data/S22/TCGA_Patients.rds")
TCGA_Patients = data.frame(TCGA_Patients)
colnames(TCGA_Patients) = c("Cancer_type", "Patient", "Drug","Response","Cancer_type_Perfect_name")

N_Cancer = readRDS("Processed_data/S22/Number_of_each_Cancer_TCGA.rds")
low_sample_drugs = c(1,3,4,5,6,7,12,13,14,19,20,21,22,23,24,25,28,29,31,33,36,37
                     ,39,40,41,42,43,45,47,49,51,52,56,57,58)
res_TCGA = readRDS("Processed_data/S24/Drug_response_TCGA_binarized.rds")
res_TCGA = res_TCGA[,-low_sample_drugs]

N_drug = ncol(res_TCGA)
drugs = data.frame(colnames(res_TCGA))

NO_tissues = list()
for (i in drugs[,1]){
  
  y = res_TCGA[!is.na(res_TCGA[,i]),i]
  length(y)
  y = data.frame(y)
  
  yy = TCGA_Patients[TCGA_Patients$Patient %in% rownames(y) & TCGA_Patients$Drug == i,]
  yy = yy[!duplicated(yy$Patient),]
  rownames(yy) = yy$Patient
  
  patient_cancer_drug_res = yy[rownames(y),]
  patient_cancer_drug_res$binarized_response = y[,1]
  patient_cancer_res = patient_cancer_drug_res[,c(1,6)]
  
  Table = table(patient_cancer_res)
  Sum = cbind(sum(Table[,1]),sum(Table[,2]))
  Table = rbind(Table, Sum)

  NO_tissues[i] = list(Table)
  
}    
    
saveRDS(NO_tissues,"Processed_data/S27/NO_tissues.rds")   
    
    
    
    