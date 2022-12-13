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


# 1) Landmark
# Boxplot  (Mean of Correlations for 50 runs for all drugs/all models) for Landmark
pdf("Figures/FS/Results_whithin/Result1/MeanCorr_Landmark_All_ML.pdf", height = 5, width = 7)
boxplot(cbind(Ridge_Landmark[,1], MLP_Landmark[,1], 
              Lasso_Landmark[,1], ENet_Landmark[,1], 
              RF_Landmark[,1]), names=NA, cex=.5, 
        main="MeanCorr", col = "#A0E0E4", ylim = c(-0.2,0.6))

dev.off()

# Boxplot  (SD of Correlations for 50 runs for all drugs/all models) for Landmark
pdf("Figures/FS/Results_whithin/Result1/STDCorr_Landmark_All_ML.pdf", height = 5, width = 7)
boxplot(cbind(Ridge_Landmark[,2], MLP_Landmark[,2], 
              Lasso_Landmark[,2], ENet_Landmark[,2], 
              RF_Landmark[,2]), names=NA, cex=.5, 
        main="STDCorr", col = "#A0E0E4")

dev.off()

# Boxplot  (Mean of MSE for 50 runs for all drugs/all models) for Landmark
pdf("Figures/FS/Results_whithin/Result1/MeanMSE_Landmark_All_ML.pdf", height = 5, width = 7)
boxplot(cbind(Ridge_Landmark[,3], MLP_Landmark[,3], 
              Lasso_Landmark[,3], ENet_Landmark[,3], 
              RF_Landmark[,3]), names=NA, cex=.5, 
        main="MeanMSE", col = "#A0E0E4")

dev.off()

# Boxplot  (SD of MSE for 50 runs for all drugs/all models) for Landmark
pdf("Figures/FS/Results_whithin/Result1/STDMSE_Landmark_All_ML.pdf", height = 5, width = 7)
boxplot(cbind(Ridge_Landmark[,4], MLP_Landmark[,4], 
              Lasso_Landmark[,4], ENet_Landmark[,4], 
              RF_Landmark[,4]), names=NA, cex=.5, 
        main="STDMSE", col = "#A0E0E4")

dev.off()



# 2) Whole Genes
# Boxplot  (Mean of Correlations for 50 runs for all drugs/all models) for Whole Genes
pdf("Figures/FS/Results_whithin/Result1/MeanCorr_WholeGenes_All_ML.pdf", height = 5, width = 7)
boxplot(cbind(Ridge_WG[,1], MLP_WG[,1], 
              Lasso_WG[,1], ENet_WG[,1], 
              RF_WG[,1]), names=NA, cex=.5, 
        main="MeanCorr", col = "#F2E6D6", ylim = c(-0.2,0.6))

dev.off()

# Boxplot  (SD of Correlations for 50 runs for all drugs/all models) for Whole Genes
pdf("Figures/FS/Results_whithin/Result1/STDCorr_WholeGenes_All_ML.pdf", height = 5, width = 7)
boxplot(cbind(Ridge_WG[,2], MLP_WG[,2], 
              Lasso_WG[,2], ENet_WG[,2], 
              RF_WG[,2]), names=NA, cex=.5, 
        main="STDCorr", col = "#F2E6D6")

dev.off()

# Boxplot  (Mean of MSE for 50 runs for all drugs/all models) for Whole Genes
pdf("Figures/FS/Results_whithin/Result1/MeanMSE_WholeGenes_All_ML.pdf", height = 5, width = 7)
boxplot(cbind(Ridge_WG[,3], MLP_WG[,3], 
              Lasso_WG[,3], ENet_WG[,3], 
              RF_WG[,3]), names=NA, cex=.5, 
        main="MeanMSE", col = "#F2E6D6")

dev.off()

# Boxplot  (SD of MSE for 50 runs for all drugs/all models) for Whole Genes
pdf("Figures/FS/Results_whithin/Result1/STDMSE_WholeGenes_All_ML.pdf", height = 5, width = 7)
boxplot(cbind(Ridge_WG[,4], MLP_WG[,4], 
              Lasso_WG[,4], ENet_WG[,4], 
              RF_WG[,4]), names=NA, cex=.5, 
        main="STDMSE", col = "#F2E6D6")

dev.off()



# 3) Drug Pathways
# Boxplot  (Mean of Correlations for 50 runs for all drugs/all models) for Drug Pathways
pdf("Figures/FS/Results_whithin/Result1/MeanCorr_DrugPathways_All_ML.pdf", height = 5, width = 7)
boxplot(cbind(Ridge_PW[,1], MLP_PW[,1], 
              Lasso_PW[,1], ENet_PW[,1], 
              RF_PW[,1]), names=NA, cex=.5, 
        main="MeanCorr", col = "#FCC0C5", ylim = c(-0.2,0.6))

dev.off()

# Boxplot  (SD of Correlations for 50 runs for all drugs/all models) for Drug Pathways
pdf("Figures/FS/Results_whithin/Result1/STDCorr_DrugPathways_All_ML.pdf", height = 5, width = 7)
boxplot(cbind(Ridge_PW[,2], MLP_PW[,2], 
              Lasso_PW[,2], ENet_PW[,2], 
              RF_PW[,2]), names=NA, cex=.5, 
        main="STDCorr", col = "#FCC0C5")

dev.off()

# Boxplot  (Mean of MSE for 50 runs for all drugs/all models) for Drug Pathways
pdf("Figures/FS/Results_whithin/Result1/MeanMSE_DrugPathways_All_ML.pdf", height = 5, width = 7)
boxplot(cbind(Ridge_PW[,3], MLP_PW[,3], 
              Lasso_PW[,3], ENet_PW[,3], 
              RF_PW[,3]), names=NA, cex=.5, 
        main="MeanMSE", col = "#FCC0C5")

dev.off()

# Boxplot  (SD of MSE for 50 runs for all drugs/all models) for Drug Pathways
pdf("Figures/FS/Results_whithin/Result1/STDMSE_DrugPathways_All_ML.pdf", height = 5, width = 7)
boxplot(cbind(Ridge_PW[,4], MLP_PW[,4], 
              Lasso_PW[,4], ENet_PW[,4], 
              RF_PW[,4]), names=NA, cex=.5, 
        main="STDMSE", col = "#FCC0C5")

dev.off()



# 4) Pathway Activities
# Boxplot  (Mean of Correlations for 50 runs for all drugs/all models) for Pathway Activities
pdf("Figures/FS/Results_whithin/Result1/MeanCorr_PathwayActivities_All_ML.pdf", height = 5, width = 7)
boxplot(cbind(Ridge_PA[,1], MLP_PA[,1], 
              Lasso_PA[,1], ENet_PA[,1], 
              RF_PA[,1]), names=NA, cex=.5, 
        main="MeanCorr", col = "#A49393", ylim = c(-0.2,0.6))

dev.off()

# Boxplot  (SD of Correlations for 50 runs for all drugs/all models) for Pathway Activities
pdf("Figures/FS/Results_whithin/Result1/STDCorr_PathwayActivities_All_ML.pdf", height = 5, width = 7)
boxplot(cbind(Ridge_PA[,2], MLP_PA[,2], 
              Lasso_PA[,2], ENet_PA[,2], 
              RF_PA[,2]), names=NA, cex=.5, 
        main="STDCorr", col = "#A49393")

dev.off()

# Boxplot  (Mean of MSE for 50 runs for all drugs/all models) for Pathway Activities
pdf("Figures/FS/Results_whithin/Result1/MeanMSE_PathwayActivities_All_ML.pdf", height = 5, width = 7)
boxplot(cbind(Ridge_PA[,3], MLP_PA[,3], 
              Lasso_PA[,3], ENet_PA[,3], 
              RF_PA[,3]), names=NA, cex=.5, 
        main="MeanMSE", col = "#A49393")

dev.off()

# Boxplot  (SD of MSE for 50 runs for all drugs/all models) for Pathway Activities
pdf("Figures/FS/Results_whithin/Result1/STDMSE_PathwayActivities_All_ML.pdf", height = 5, width = 7)
boxplot(cbind(Ridge_PA[,4], MLP_PA[,4], 
              Lasso_PA[,4], ENet_PA[,4], 
              RF_PA[,4]), names=NA, cex=.5, 
        main="STDMSE", col = "#A49393")

dev.off()



# 5) Transcription Factors
# Boxplot  (Mean of Correlations for 50 runs for all drugs/all models) for Transcription Factors
pdf("Figures/FS/Results_whithin/Result1/MeanCorr_TranscriptionFactors_All_ML.pdf", height = 5, width = 7)
  boxplot(cbind(Ridge_TF[,1], MLP_TF[,1], 
              Lasso_TF[,1], ENet_TF[,1], 
              RF_TF[,1]), names=NA, cex=.5, 
        main="MeanCorr", col = "#F582A8", ylim = c(-0.2,0.6))

dev.off()

# Boxplot  (SD of Correlations for 50 runs for all drugs/all models) for Transcription Factors
pdf("Figures/FS/Results_whithin/Result1/STDCorr_TranscriptionFactors_All_ML.pdf", height = 5, width = 7)
boxplot(cbind(Ridge_TF[,2], MLP_TF[,2], 
              Lasso_TF[,2], ENet_TF[,2], 
              RF_TF[,2]), names=NA, cex=.5, 
        main="STDCorr", col = "#F582A8")

dev.off()

# Boxplot  (Mean of MSE for 50 runs for all drugs/all models) for Transcription Factors
pdf("Figures/FS/Results_whithin/Result1/MeanMSE_TranscriptionFactors_All_ML.pdf", height = 5, width = 7)
boxplot(cbind(Ridge_TF[,3], MLP_TF[,3], 
              Lasso_TF[,3], ENet_TF[,3], 
              RF_TF[,3]), names=NA, cex=.5, 
        main="MeanMSE", col = "#F582A8")

dev.off()

# Boxplot  (SD of MSE for 50 runs for all drugs/all models) for Transcription Factors
pdf("Figures/FS/Results_whithin/Result1/STDMSE_TranscriptionFactors_All_ML.pdf", height = 5, width = 7)
boxplot(cbind(Ridge_TF[,4], MLP_TF[,4], 
              Lasso_TF[,4], ENet_TF[,4], 
              RF_TF[,4]), names=NA, cex=.5, 
        main="STDMSE", col = "#F582A8")

dev.off()


