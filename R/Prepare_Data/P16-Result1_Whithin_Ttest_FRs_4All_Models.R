rm(list=ls())

library(ggplot2)
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
sen = readRDS("Processed_data/S1/sensitivity_matrix_AUC.rds")
N_drugs = 1448

#Read Landmark results
RF_Landmark = c()
ENet_Landmark = c()
Lasso_Landmark = c()
Ridge_Landmark = c()
MLP_Landmark = c()

for(i in 1:N_drugs){
  print(i)
  R_Landmark = readRDS(paste0("Processed_from_SLURM/Results_Landmark_All_Models/Result_",as.character(i),".rds"))
  Ridge_Landmark = rbind(Ridge_Landmark, R_Landmark[4,])
  MLP_Landmark = rbind(MLP_Landmark, R_Landmark[5,])
  Lasso_Landmark = rbind(Lasso_Landmark, R_Landmark[3,])
  ENet_Landmark = rbind(ENet_Landmark, R_Landmark[2,])
  RF_Landmark = rbind(RF_Landmark, R_Landmark[1,])
}

rownames(Ridge_Landmark) = colnames(sen)
rownames(MLP_Landmark) = colnames(sen)
rownames(Lasso_Landmark) = colnames(sen)
rownames(ENet_Landmark) = colnames(sen)
rownames(RF_Landmark) = colnames(sen)


#Read Whole Genes results
RF_WG = c()
ENet_WG = c()
Lasso_WG = c()
Ridge_WG = c()
MLP_WG = c()

for(i in 1:N_drugs){
  print(i)
  R_WG = readRDS(paste0("Processed_from_SLURM/Results_WholeGenes_All_Models/Result_",as.character(i),".rds"))
  Ridge_WG = rbind(Ridge_WG, R_WG[4,])
  MLP_WG = rbind(MLP_WG, R_WG[5,])
  Lasso_WG = rbind(Lasso_WG, R_WG[3,])
  ENet_WG = rbind(ENet_WG, R_WG[2,])
  RF_WG = rbind(RF_WG, R_WG[1,])
}

rownames(Ridge_WG) = colnames(sen)
rownames(MLP_WG) = colnames(sen)
rownames(Lasso_WG) = colnames(sen)
rownames(ENet_WG) = colnames(sen)
rownames(RF_WG) = colnames(sen)



#Read Drug Pathway results
RF_PW = c()
ENet_PW = c()
Lasso_PW = c()
Ridge_PW = c()
MLP_PW = c()
I_zeros = c()
for(i in 1:N_drugs){
  print(i)
  R_PW = readRDS(paste0("Processed_from_SLURM/Results_DrugPathways_All_Models/Result_",as.character(i),".rds"))
  if(is.null(nrow(R_PW))){
    I_zeros = c(I_zeros,i)
  }else{
    Ridge_PW = rbind(Ridge_PW, R_PW[4,])
    MLP_PW = rbind(MLP_PW, R_PW[5,])
    Lasso_PW = rbind(Lasso_PW, R_PW[3,])
    ENet_PW = rbind(ENet_PW, R_PW[2,])
    RF_PW = rbind(RF_PW, R_PW[1,])
  }
}

Ridge_PW_new = c()
MLP_PW_new = c()
Lasso_PW_new = c()
ENet_PW_new = c()
RF_PW_new = c()

C = 1
for(l in 1:N_drugs){
  if(l %in% I_zeros){
    
    Ridge_PW_new = rbind(Ridge_PW_new, rep(0,5))
    MLP_PW_new = rbind(MLP_PW_new, rep(0,5))
    Lasso_PW_new = rbind(Lasso_PW_new, rep(0,5))
    ENet_PW_new = rbind(ENet_PW_new, rep(0,5))
    RF_PW_new = rbind(RF_PW_new, rep(0,5))
    
  }else{
    
    Ridge_PW_new = rbind(Ridge_PW_new, Ridge_PW[C,])
    MLP_PW_new = rbind(MLP_PW_new, MLP_PW[C,])
    Lasso_PW_new = rbind(Lasso_PW_new, Lasso_PW[C,])
    ENet_PW_new = rbind(ENet_PW_new, ENet_PW[C,])
    RF_PW_new = rbind(RF_PW_new, RF_PW[C,])
    
    C = C+1
  }
}

rownames(Ridge_PW_new) = colnames(sen)
rownames(MLP_PW_new) = colnames(sen)
rownames(Lasso_PW_new) = colnames(sen)
rownames(ENet_PW_new) = colnames(sen)
rownames(RF_PW_new) = colnames(sen)


#Read Progeny Pathway activities results
RF_PA = c()
ENet_PA = c()
Lasso_PA = c()
Ridge_PA = c()
MLP_PA = c()
I = c()
for(i in 1:N_drugs){
  print(i)
  if(file.exists(paste0("Processed_from_SLURM/Results_Progeny_Pathway_Activities_All_Models/Result_",as.character(i),".rds"))){
    R_PA = readRDS(paste0("Processed_from_SLURM/Results_Progeny_Pathway_Activities_All_Models/Result_",as.character(i),".rds"))
    Ridge_PA = rbind(Ridge_PA, R_PA[4,])
    MLP_PA = rbind(MLP_PA, R_PA[5,])
    Lasso_PA = rbind(Lasso_PA, R_PA[3,])
    ENet_PA = rbind(ENet_PA, R_PA[2,])
    RF_PA = rbind(RF_PA, R_PA[1,])
  }else{
    I = c(I,i)
  }
}

Ridge_PA_new = c()
MLP_PA_new = c()
Lasso_PA_new = c()
ENet_PA_new = c()
RF_PA_new = c()

C = 1
for(l in 1:N_drugs){
  if(l %in% I){
    Ridge_PA_new = rbind(Ridge_PA_new, rep(0,4))
    MLP_PA_new = rbind(MLP_PA_new, rep(0,4))
    Lasso_PA_new = rbind(Lasso_PA_new, rep(0,4))
    ENet_PA_new = rbind(ENet_PA_new, rep(0,4))
    RF_PA_new = rbind(RF_PA_new, rep(0,4))
    
  }else{
    Ridge_PA_new = rbind(Ridge_PA_new, Ridge_PA[C,])
    MLP_PA_new = rbind(MLP_PA_new, MLP_PA[C,])
    Lasso_PA_new = rbind(Lasso_PA_new, Lasso_PA[C,])
    ENet_PA_new = rbind(ENet_PA_new, ENet_PA[C,])
    RF_PA_new = rbind(RF_PA_new, RF_PA[C,])
    
    C = C+1
  }
}

rownames(Ridge_PA_new) = colnames(sen)
rownames(MLP_PA_new) = colnames(sen)
rownames(Lasso_PA_new) = colnames(sen)
rownames(ENet_PA_new) = colnames(sen)
rownames(RF_PA_new) = colnames(sen)


#Read decoupleR Transcription factor activities results
RF_TF = c()
ENet_TF = c()
Lasso_TF = c()
Ridge_TF = c()
MLP_TF = c()

for(i in 1:N_drugs){
  print(i)
  R_TF = readRDS(paste0("Processed_from_SLURM/Results_decoupleR_TF_Activities_All_Models/Result_",as.character(i),".rds"))
  Ridge_TF = rbind(Ridge_TF, R_TF[4,])
  MLP_TF = rbind(MLP_TF, R_TF[5,])
  Lasso_TF = rbind(Lasso_TF, R_TF[3,])
  ENet_TF = rbind(ENet_TF, R_TF[2,])
  RF_TF = rbind(RF_TF, R_TF[1,])
}

rownames(Ridge_TF) = colnames(sen)
rownames(MLP_TF) = colnames(sen)
rownames(Lasso_TF) = colnames(sen)
rownames(ENet_TF) = colnames(sen)
rownames(RF_TF) = colnames(sen)


# Statistical Test FRs
## Ridge
wilcox.test(Ridge_Landmark[,1], Ridge_WG[,1], paired = FALSE, alternative = "greater")
wilcox.test(Ridge_Landmark[,1], Ridge_PW[,1], paired = FALSE, alternative = "greater")
wilcox.test(Ridge_Landmark[,1], Ridge_PA[,1], paired = FALSE, alternative = "greater")
wilcox.test(Ridge_Landmark[,1], Ridge_TF[,1], paired = FALSE, alternative = "greater")

wilcox.test(Ridge_PW[,1], Ridge_WG[,1], paired = FALSE, alternative = "greater")
wilcox.test(Ridge_PW[,1], Ridge_PA[,1], paired = FALSE, alternative = "greater")
wilcox.test(Ridge_PW[,1], Ridge_TF[,1], paired = FALSE, alternative = "greater")

wilcox.test(Ridge_WG[,1], Ridge_TF[,1], paired = FALSE, alternative = "greater")
wilcox.test(Ridge_WG[,1], Ridge_PA[,1], paired = FALSE, alternative = "greater")

wilcox.test(Ridge_TF[,1], Ridge_PA[,1], paired = FALSE, alternative = "greater")


## MLP
wilcox.test(MLP_Landmark[,1], MLP_WG[,1], paired = FALSE, alternative = "greater")
wilcox.test(MLP_Landmark[,1], MLP_PW[,1], paired = FALSE, alternative = "greater")
wilcox.test(MLP_Landmark[,1], MLP_PA[,1], paired = FALSE, alternative = "greater")
wilcox.test(MLP_Landmark[,1], MLP_TF[,1], paired = FALSE, alternative = "greater")

wilcox.test(MLP_PW[,1], MLP_WG[,1], paired = FALSE, alternative = "greater")
wilcox.test(MLP_PW[,1], MLP_PA[,1], paired = FALSE, alternative = "greater")
wilcox.test(MLP_PW[,1], MLP_TF[,1], paired = FALSE, alternative = "greater")

wilcox.test(MLP_WG[,1], MLP_TF[,1], paired = FALSE, alternative = "greater")
wilcox.test(MLP_WG[,1], MLP_PA[,1], paired = FALSE, alternative = "greater")

wilcox.test(MLP_TF[,1], MLP_PA[,1], paired = FALSE, alternative = "greater")


## RF
wilcox.test(RF_Landmark[,1], RF_WG[,1], paired = FALSE, alternative = "less")
wilcox.test(RF_Landmark[,1], RF_PW[,1], paired = FALSE, alternative = "greater")
wilcox.test(RF_Landmark[,1], RF_PA[,1], paired = FALSE, alternative = "greater")
wilcox.test(RF_Landmark[,1], RF_TF[,1], paired = TRUE, alternative = "greater")

wilcox.test(RF_PW[,1], RF_WG[,1], paired = FALSE, alternative = "greater")
wilcox.test(RF_PW[,1], RF_PA[,1], paired = FALSE, alternative = "greater")
wilcox.test(RF_PW[,1], RF_TF[,1], paired = FALSE, alternative = "greater")

wilcox.test(RF_WG[,1], RF_TF[,1], paired = FALSE, alternative = "greater")
wilcox.test(RF_WG[,1], RF_PA[,1], paired = FALSE, alternative = "greater")

wilcox.test(RF_TF[,1], RF_PA[,1], paired = FALSE, alternative = "greater")


## Lasso
wilcox.test(Lasso_Landmark[,1], Lasso_WG[,1], paired = FALSE, alternative = "less")
wilcox.test(Lasso_Landmark[,1], Lasso_PW[,1], paired = FALSE, alternative = "less")
wilcox.test(Lasso_Landmark[,1], Lasso_PA[,1], paired = FALSE, alternative = "less")
wilcox.test(Lasso_Landmark[,1], Lasso_TF[,1], paired = TRUE, alternative = "greater")

wilcox.test(Lasso_PW[,1], Lasso_WG[,1], paired = FALSE, alternative = "less")
wilcox.test(Lasso_PW[,1], Lasso_PA[,1], paired = FALSE, alternative = "greater")
wilcox.test(Lasso_PW[,1], Lasso_TF[,1], paired = FALSE, alternative = "greater")

wilcox.test(Lasso_WG[,1], Lasso_TF[,1], paired = FALSE, alternative = "greater")
wilcox.test(Lasso_WG[,1], Lasso_PA[,1], paired = FALSE, alternative = "greater")

wilcox.test(Lasso_TF[,1], Lasso_PA[,1], paired = FALSE, alternative = "less")


## ENet
wilcox.test(ENet_Landmark[,1], ENet_WG[,1], paired = FALSE, alternative = "less")
wilcox.test(ENet_Landmark[,1], ENet_PW[,1], paired = FALSE, alternative = "less")
wilcox.test(ENet_Landmark[,1], ENet_PA[,1], paired = FALSE, alternative = "greater")
wilcox.test(ENet_Landmark[,1], ENet_TF[,1], paired = TRUE, alternative = "greater")

wilcox.test(ENet_PW[,1], ENet_WG[,1], paired = FALSE, alternative = "less")
wilcox.test(ENet_PW[,1], ENet_PA[,1], paired = FALSE, alternative = "greater")
wilcox.test(ENet_PW[,1], ENet_TF[,1], paired = FALSE, alternative = "greater")

wilcox.test(ENet_WG[,1], ENet_TF[,1], paired = FALSE, alternative = "greater")
wilcox.test(ENet_WG[,1], ENet_PA[,1], paired = FALSE, alternative = "greater")

wilcox.test(ENet_TF[,1], ENet_PA[,1], paired = FALSE, alternative = "less")



# Statistical Test Models
# Landmark
t.test(Ridge_Landmark[,1], MLP_Landmark[,1], paired = FALSE, alternative = "greater")
t.test(Ridge_Landmark[,1], RF_Landmark[,1], paired = FALSE, alternative = "greater")
t.test(Ridge_Landmark[,1], Lasso_Landmark[,1], paired = FALSE, alternative = "greater")
t.test(Ridge_Landmark[,1], ENet_Landmark[,1], paired = FALSE, alternative = "greater")

t.test(MLP_Landmark[,1], RF_Landmark[,1], paired = FALSE, alternative = "greater")
t.test(MLP_Landmark[,1], Lasso_Landmark[,1], paired = FALSE, alternative = "greater")
t.test(MLP_Landmark[,1], ENet_Landmark[,1], paired = FALSE, alternative = "greater")

t.test(RF_Landmark[,1], Lasso_Landmark[,1], paired = FALSE, alternative = "greater")
t.test(RF_Landmark[,1], ENet_Landmark[,1], paired = FALSE, alternative = "greater")

t.test(ENet_Landmark[,1], Lasso_Landmark[,1], paired = FALSE, alternative = "greater")




# Whole genes
t.test(Ridge_WG[,1], MLP_WG[,1], paired = FALSE, alternative = "greater")
t.test(Ridge_WG[,1], RF_WG[,1], paired = FALSE, alternative = "greater")
t.test(Ridge_WG[,1], Lasso_WG[,1], paired = FALSE, alternative = "greater")
t.test(Ridge_WG[,1], ENet_WG[,1], paired = FALSE, alternative = "greater")


t.test(MLP_WG[,1], RF_WG[,1], paired = FALSE, alternative = "less")
t.test(MLP_WG[,1], Lasso_WG[,1], paired = FALSE, alternative = "greater")
t.test(MLP_WG[,1], ENet_WG[,1], paired = FALSE, alternative = "greater")

t.test(RF_WG[,1], Lasso_WG[,1], paired = FALSE, alternative = "greater")
t.test(RF_WG[,1], ENet_WG[,1], paired = FALSE, alternative = "greater")

t.test(ENet_WG[,1], Lasso_WG[,1], paired = FALSE, alternative = "greater")


# Drug Pathway
t.test(Ridge_PW[,1], MLP_PW[,1], paired = FALSE, alternative = "greater")
t.test(Ridge_PW[,1], RF_PW[,1], paired = FALSE, alternative = "greater")
t.test(Ridge_PW[,1], Lasso_PW[,1], paired = FALSE, alternative = "greater")
t.test(Ridge_PW[,1], ENet_PW[,1], paired = FALSE, alternative = "greater")

t.test(MLP_PW[,1], RF_PW[,1], paired = FALSE, alternative = "greater")
t.test(MLP_PW[,1], Lasso_PW[,1], paired = FALSE, alternative = "greater")
t.test(MLP_PW[,1], ENet_PW[,1], paired = FALSE, alternative = "greater")

t.test(RF_PW[,1], Lasso_PW[,1], paired = FALSE, alternative = "greater")
t.test(RF_PW[,1], ENet_PW[,1], paired = FALSE, alternative = "greater")

t.test(ENet_PW[,1], Lasso_PW[,1], paired = FALSE, alternative = "greater")

#Pathway activities
t.test(Ridge_PA[,1], MLP_PA[,1], paired = FALSE, alternative = "greater")
t.test(Ridge_PA[,1], RF_PA[,1], paired = FALSE, alternative = "greater")
t.test(Ridge_PA[,1], Lasso_PA[,1], paired = FALSE, alternative = "greater")
t.test(Ridge_PA[,1], ENet_PA[,1], paired = FALSE, alternative = "greater")

t.test(MLP_PA[,1], RF_PA[,1], paired = FALSE, alternative = "greater")
t.test(MLP_PA[,1], Lasso_PA[,1], paired = FALSE, alternative = "greater")
t.test(MLP_PA[,1], ENet_PA[,1], paired = FALSE, alternative = "greater")

t.test(RF_PA[,1], Lasso_PA[,1], paired = FALSE, alternative = "greater")
t.test(RF_PA[,1], ENet_PA[,1], paired = FALSE, alternative = "greater")

t.test(ENet_PA[,1], Lasso_PA[,1], paired = FALSE, alternative = "greater")


# Transcription Factors
t.test(Ridge_TF[,1], MLP_TF[,1], paired = FALSE, alternative = "greater")
t.test(Ridge_TF[,1], RF_TF[,1], paired = FALSE, alternative = "greater")
t.test(Ridge_TF[,1], Lasso_TF[,1], paired = FALSE, alternative = "greater")
t.test(Ridge_TF[,1], ENet_TF[,1], paired = FALSE, alternative = "greater")

t.test(MLP_TF[,1], RF_TF[,1], paired = FALSE, alternative = "greater")
t.test(MLP_TF[,1], Lasso_TF[,1], paired = FALSE, alternative = "greater")
t.test(MLP_TF[,1], ENet_TF[,1], paired = FALSE, alternative = "greater")

t.test(RF_TF[,1], Lasso_TF[,1], paired = FALSE, alternative = "greater")
t.test(RF_TF[,1], ENet_TF[,1], paired = FALSE, alternative = "greater")

t.test(ENet_TF[,1], Lasso_TF[,1], paired = FALSE, alternative = "less")







