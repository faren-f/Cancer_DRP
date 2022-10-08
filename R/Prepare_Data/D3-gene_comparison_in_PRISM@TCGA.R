rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
sen_PRISM = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")
res_TCGA = readRDS("Processed_data/Other/Res_TCGA_24_Drugs.rds")

GE = readRDS("Processed_data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")
GE_TCGA = readRDS("Processed_data/S23/expresion_matrix_TCGA.rds")


# Remove genes whose Q3 is zero
q3_genes = apply(GE_TCGA,2,quantile,prob=0.75)
GE_TCGA = GE_TCGA[,-which(q3_genes==0)]
GE = GE[,-which(q3_genes==0)]

N_drug = ncol(sen_PRISM)
drugs = data.frame(colnames(sen_PRISM))
#saveRDS(drugs,"Processed_data/Other/24_drugs.rds")
d = 23
drug = drugs[d,1]


#source("F21-Drug_Pathway_Level_genes.R")
#i=1
#pathway_gene_set = Drug_Pathway_gene_set(drug = drug, level=i)

# I = intersect(colnames(GE),pathway_gene_set)
# X = GE[,I]
# X_TCGA = GE_TCGA[,I]
X = GE
X_TCGA = GE_TCGA
source("F15-Feature_Selection_PRISM@TCGA.R")
selected_features = c("Landmark_genes")
Omics_List = Feature_Selection_PRISM_TCGA(selected_features,GE = X ,GE_test = X_TCGA)
X = Omics_List[[1]]
X_TCGA = Omics_List[[3]]


X_P = X[!is.na(sen_PRISM[,d]),]
y_P = sen_PRISM[!is.na(sen_PRISM[,d]),d]

X_T = X_TCGA[!is.na(res_TCGA[,d]),]
y_T = res_TCGA[!is.na(res_TCGA[,d]),d]

TCGA_Patients = readRDS("Processed_data/S22/TCGA_Patients.rds")
s = TCGA_Patients[TCGA_Patients[,3]=="carboplatin",]
#table(s[,5])
#Ti =s[s[,1]=="UCS"|s[,1]== "UCEC",2]
Ti =s[!(s[,1]=="HNSC"),2]


y_T = y_T[rownames(X_T)%in%Ti]
X_T = X_T[rownames(X_T)%in%Ti,]

source("F18-Combat_Normalization.R")
X_Normalization = Combat_Scale(X_P,X_T)

X_P = X_Normalization[[1]]
X_T = X_Normalization[[2]]
N_genes = ncol(X_P)

hist(y_P)

# Cor_P = abs(cor(X_P,y_P))
# sort(Cor_P, decreasing = TRUE)[1:50]
# order_P = order(Cor_P, decreasing = TRUE)[1:50]
# 
# 
# Cor_T = abs(cor(X_T,y_T))
# sort(Cor_T, decreasing = TRUE)[1:50]
# order_T = order(Cor_T, decreasing = TRUE)[1:50]
# 
# length(intersect(order_P,order_T))



X_P_nonres = X_P[y_P>1.2,]
X_P_res = X_P[y_P<1,]
X_T = t(X_T)

X_P_nonres = t(X_P_nonres)
X_T_nonres = cor(X_T,X_P_nonres)

X_P_res = t(X_P_res)
X_T_res = cor(X_T,X_P_res)

X_T_nonres_ave = apply(X_T_nonres,1,mean)
X_T_res_ave = apply(X_T_res,1,mean)

Respond = rep(0,length(X_T_nonres_ave))
for(i in 1:length(X_T_nonres_ave)){
  Respond[i] = ifelse(X_T_nonres_ave[i] > X_T_res_ave[i], 1, 2)
}

C = data.frame(cbind(Respond,y_T))
table(C)

chisq.test(table(C), simulate.p.value = TRUE)




