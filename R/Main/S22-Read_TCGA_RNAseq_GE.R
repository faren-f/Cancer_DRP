rm(list=ls())
library("readxl")

# Read clinical data from table (2016-Ding)
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
res_TCGA = readRDS("Processed_data/S21/Drug_Response_matrix_TCGA.rds")
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

Cancer_type = rep(Cancer,times = S)
Response_table$Cancer_type = Cancer_type
TCGA_Patients = cbind(Cancer_type,Response_table$bcr_patient_barcode,
                    Response_table$drug_name, Response_table$measure_of_response,
                    Response_table$Cancer)
## Read RNAseq data from (http://firebrowse.org) 

ACC = read.table("Raw_data/TCGA/RNAseq/ACC/ACC.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
ACC = ACC[,1:(ncol(ACC)-1)]
BLCA = read.table("Raw_data/TCGA/RNAseq/BLCA/BLCA.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
BRCA = read.table("Raw_data/TCGA/RNAseq/BRCA/BRCA.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
BRCA = BRCA[,1:(ncol(BRCA)-1)]
CESC = read.table("Raw_data/TCGA/RNAseq/CESC/CESC.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
#CHOL = read.table("Raw_data/TCGA/RNAseq/CHOL/CHOL.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
COAD = read.table("Raw_data/TCGA/RNAseq/COAD/COAD.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
#COADREAD = read.table("Raw_data/TCGA/RNAseq/COADREAD/COADREAD.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
#DLBC = read.table("Raw_data/TCGA/RNAseq/DLBC/DLBC.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
ESCA = read.table("Raw_data/TCGA/RNAseq/ESCA/ESCA.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
#GBM = read.table("Raw_data/TCGA/RNAseq/GBM/GBM.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
#GBMLGG = read.table("Raw_data/TCGA/RNAseq/GBMLGG/GBMLGG.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
HNSC = read.table("Raw_data/TCGA/RNAseq/HNSC/HNSC.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
#KICH = read.table("Raw_data/TCGA/RNAseq/KICH/KICH.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
#KIPAN = read.table("Raw_data/TCGA/RNAseq/KIPAN/KIPAN.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
KIRC = read.table("Raw_data/TCGA/RNAseq/KIRC/KIRC.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
KIRP = read.table("Raw_data/TCGA/RNAseq/KIRP/KIRP.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
#LAML = read.table("Raw_data/TCGA/RNAseq/LAML/LAML.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
LGG = read.table("Raw_data/TCGA/RNAseq/LGG/LGG.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
LIHC = read.table("Raw_data/TCGA/RNAseq/LIHC/LIHC.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
LUAD = read.table("Raw_data/TCGA/RNAseq/LUAD/LUAD.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
LUSC = read.table("Raw_data/TCGA/RNAseq/LUSC/LUSC.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
MESO = read.table("Raw_data/TCGA/RNAseq/MESO/MESO.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
OV = read.table("Raw_data/TCGA/RNAseq/OV/OV.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
PAAD = read.table("Raw_data/TCGA/RNAseq/PAAD/PAAD.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
PAAD = PAAD[,1:(ncol(PAAD)-1)]
PCPG = read.table("Raw_data/TCGA/RNAseq/PCPG/PCPG.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",fill = TRUE, header=FALSE)
PRAD = read.table("Raw_data/TCGA/RNAseq/PRAD/PRAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",fill = TRUE, header=FALSE)
READ = read.table("Raw_data/TCGA/RNAseq/READ/READ.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
SARC = read.table("Raw_data/TCGA/RNAseq/SARC/SARC.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
SKCM = read.table("Raw_data/TCGA/RNAseq/SKCM/SKCM.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
STAD = read.table("Raw_data/TCGA/RNAseq/STAD/STAD.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
#STES = read.table("Raw_data/TCGA/RNAseq/STES/STES.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
TGCT = read.table("Raw_data/TCGA/RNAseq/TGCT/TGCT.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
THCA = read.table("Raw_data/TCGA/RNAseq/THCA/THCA.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
#THYM = read.table("Raw_data/TCGA/RNAseq/THYM/THYM.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
UCEC = read.table("Raw_data/TCGA/RNAseq/UCEC/UCEC.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
UCS = read.table("Raw_data/TCGA/RNAseq/UCS/UCS.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
#UVM = read.table("Raw_data/TCGA/RNAseq/UVM/UVM.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)

Cancers = c("ACC","BLCA","BRCA","CESC","COAD","ESCA",
            "HNSC","KIRC","KIRP","LGG","LIHC","LUAD",
            "LUSC","MESO","OV","PAAD","PCPG","PRAD","READ",
            "SARC","SKCM","STAD","TGCT","THCA","UCEC","UCS")

TCGA_GE = c()
TT = c()
for (i in Cancers){
  print(i)
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
  
  colnames(Cancer_i) = Cancer_i[1,]
  Cancer_i = Cancer_i[-1,]
  rownames(Cancer_i) = Cancer_i[,1]
  Cancer_i = Cancer_i[,c(-1)]
  Cancer_i = t(Cancer_i)
  #B = data.frame((rownames(Cancer_i)))
  Cancer_i_2 = apply(Cancer_i,2,as.numeric)
  
  if(length(dup)>0){
    dup = dup-1
    for (k in dup){
      rep = apply(Cancer_i_2[c(k-1,k),],2,mean)
      Cancer_i_2[k-1,] = rep
    }
    Cancer_i_2 = Cancer_i_2[-dup,]
    Cancer_i = Cancer_i[-dup,]
  }
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

N_TT = c()
for(u in unique(TT)){
  N_TT = c(N_TT, sum(TT %in% u))
}

Number_of_each_Cancer = cbind(unique(TT),N_TT)

I_samples = intersect(rownames(TCGA_GE),rownames(res_TCGA))

TCGA_GE = TCGA_GE[I_samples,]
res_TCGA = res_TCGA[I_samples,]

saveRDS(Number_of_each_Cancer,"Processed_data/S22/Number_of_each_Cancer_TCGA.rds")
saveRDS(TCGA_GE,"Processed_data/S22/expresion_matrix_TCGA.rds")
saveRDS(res_TCGA,"Processed_data/S22/Drug_response_matrix_TCGA.rds")
saveRDS(TCGA_Patients,"Processed_data/S22/TCGA_Patients.rds")
