rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
TCGA_Patients = readRDS("Processed_data/S22/TCGA_Patients.rds")
TCGA_Patients = data.frame(TCGA_Patients)
colnames(TCGA_Patients) = c("Cancer_type", "Patient", "Drug","Response","Cancer_type_Perfect_name")

N_Cancer = readRDS("Processed_data/S22/Number_of_each_Cancer_TCGA.rds")

res_TCGA = readRDS("Processed_data/Other/Res_TCGA_24_Drugs.rds")
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
    
    
    
    