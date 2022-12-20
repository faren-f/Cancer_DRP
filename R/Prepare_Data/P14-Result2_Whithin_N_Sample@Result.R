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


# Finding Number of samples for each drug
N_Na = data.frame(apply(sen,2,function(x){sum(is.na(x))}))
N_S = data.frame(475-N_Na[,1])
rownames(N_S) = rownames(N_Na)

#1) Landmark
# Plot number of samples against  results for all the ML methods

#1-1) Ridge
df_Ridge_Landmark = data.frame(Cor = Ridge_Landmark[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/Result2/No_Samples&Cor_Landmark_Ridge.pdf", height = 5, width = 6)

ggplot(df_Ridge_Landmark, aes(x = Cor, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE)+ theme_classic()

dev.off()
lm(Cor~N, data=df_Ridge_Landmark)
cor(df_Ridge_Landmark$Cor,df_Ridge_Landmark$N)

#1-2) MLP
df_MLP_Landmark = data.frame(Cor = MLP_Landmark[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/Result2/No_Samples&Cor_Landmark_MLP.pdf",height = 5, width = 6)

ggplot(df_MLP_Landmark, aes(x = Cor, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) + theme_classic()

dev.off()

#1-3) Lasso
df_Lasso_Landmark = data.frame(Cor = Lasso_Landmark[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/Result2/No_Samples&Cor_Landmark_Lasso.pdf",height = 5, width = 6)

ggplot(df_Lasso_Landmark, aes(x = Cor, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) + theme_classic()

dev.off()

#1-4) ENet
df_ENet_Landmark = data.frame(Cor = ENet_Landmark[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/Result2/No_Samples&Cor_Landmark_ENet.pdf",height = 5, width = 6)
ggplot(df_ENet_Landmark, aes(x = Cor, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) + theme_classic()

dev.off()

#1-5) RF
df_RF_Landmark = data.frame(Cor = RF_Landmark[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/Result2/No_Samples&Cor_Landmark_RF.pdf",height = 5, width = 6)
ggplot(df_RF_Landmark, aes(x = Cor, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) + theme_classic()

dev.off()



#2) WholeGenes
# Plot number of samples against  results for all the ML methods

#1-1) Ridge
df_Ridge_WG = data.frame(Cor = Ridge_WG[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/Result2/No_Samples&Cor_WholeGenes_Ridge.pdf",height = 5, width = 6)

ggplot(df_Ridge_WG, aes(x = Cor, y = N)) + geom_point(size = 0.5) + 
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) + theme_classic()

dev.off()
cor(df_Ridge_WG$Cor,df_Ridge_WG$N)

lm(Cor~N, data=df_Ridge_WG)

#1-2) MLP
df_MLP_WG = data.frame(Cor = MLP_WG[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/Result2/No_Samples&Cor_WholeGenes_MLP.pdf",height = 5, width = 6)

ggplot(df_MLP_WG, aes(x = Cor, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) + theme_classic()

dev.off()

#1-3) Lasso
df_Lasso_WG = data.frame(Cor = Lasso_WG[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/Result2/No_Samples&Cor_WholeGenes_Lasso.pdf",height = 5, width = 6)

ggplot(df_Lasso_WG, aes(x = Cor, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) + theme_classic()

dev.off()

#1-4) ENet
df_ENet_WG = data.frame(Cor = ENet_WG[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/Result2/No_Samples&Cor_WholeGenes_ENet.pdf",height = 5, width = 6)
ggplot(df_ENet_WG, aes(x = Cor, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) + theme_classic()

dev.off()

#1-5) RF
df_RF_WG = data.frame(Cor = RF_WG[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/Result2/No_Samples&Cor_WholeGenes_RF.pdf",height = 5, width = 6)
ggplot(df_RF_WG, aes(x = Cor, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) + theme_classic()

dev.off()



#3) Drug Pathways
# Plot number of samples against  results for all the ML methods

#1-1) Ridge
df_Ridge_PW = data.frame(Cor = Ridge_PW_new[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/Result2/No_Samples&Cor_DrugPathways_Ridge.pdf",height = 5, width = 6)

ggplot(df_Ridge_PW, aes(x = Cor, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE)+ theme_classic()

dev.off()

R1 = Ridge_PW[Ridge_PW$N_genes<300,]
cor(R1$Mean_Corr,R1$N_genes)

R2 = Ridge_PW[Ridge_PW$N_genes<500,]
cor(R2$Mean_Corr,R2$N_genes)

R3 = Ridge_PW[Ridge_PW$N_genes<700,]
cor(R3$Mean_Corr,R3$N_genes)

R4 = Ridge_PW[Ridge_PW$N_genes<1000,]
cor(R4$Mean_Corr,R4$N_genes)

R5 = Ridge_PW[Ridge_PW$N_genes<1500,]
cor(R5$Mean_Corr,R5$N_genes)

R6 = Ridge_PW[Ridge_PW$N_genes>1,]
cor(R6$Mean_Corr,R6$N_genes)


cor(Ridge_PW$Mean_Corr,Ridge_PW$N_genes)
lm(Mean_Corr~N_genes,Ridge_PW)
#1-2) MLP
df_MLP_PW = data.frame(Cor = MLP_PW[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/Result2/No_Samples&Cor_DrugPathways_MLP.pdf",height = 5, width = 6)

ggplot(df_MLP_PW, aes(x = Cor, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) + theme_classic()

dev.off()

#1-3) Lasso
df_Lasso_PW = data.frame(Cor = Lasso_PW[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/Result2/No_Samples&Cor_DrugPathways_Lasso.pdf",height = 5, width = 6)

ggplot(df_Lasso_PW, aes(x = Cor, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) + theme_classic()

dev.off()

#1-4) ENet
df_ENet_PW = data.frame(Cor = ENet_PW[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/Result2/No_Samples&Cor_DrugPathways_ENet.pdf",height = 5, width = 6)
ggplot(df_ENet_PW, aes(x = Cor, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) + theme_classic()

dev.off()

#1-5) RF
df_RF_PW = data.frame(Cor = RF_PW[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/Result2/No_Samples&Cor_DrugPathways_RF.pdf",height = 5, width = 6)
ggplot(df_RF_PW, aes(x = Cor, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) + theme_classic()

dev.off()


#4) Pathway Activities
# Plot number of samples against  results for all the ML methods

#1-1) Ridge
df_Ridge_PA = data.frame(Cor = Ridge_PA_new[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/Result2/No_Samples&Cor_PathwayActivities_Ridge.pdf",height = 5, width = 6)

ggplot(df_Ridge_PA, aes(x = Cor, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE)+ theme_classic()

dev.off()
cor(df_Ridge_PA$Cor,df_Ridge_PA$N)

lm (df_Ridge_PA)
#1-2) MLP
df_MLP_PA = data.frame(Cor = MLP_PA[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/Result2/No_Samples&Cor_PathwayActivities_MLP.pdf",height = 5, width = 6)

ggplot(df_MLP_PA, aes(x = Cor, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) + theme_classic()

dev.off()

#1-3) Lasso
df_Lasso_PA = data.frame(Cor = Lasso_PA[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/Result2/No_Samples&Cor_PathwayActivities_Lasso.pdf",height = 5, width = 6)

ggplot(df_Lasso_PA, aes(x = Cor, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) + theme_classic()

dev.off()

#1-4) ENet
df_ENet_PA = data.frame(Cor = ENet_PA[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/Result2/No_Samples&Cor_PathwayActivities_ENet.pdf",height = 5, width = 6)
ggplot(df_ENet_PA, aes(x = Cor, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) + theme_classic()

dev.off()

#1-5) RF
df_RF_PA = data.frame(Cor = RF_PA[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/Result2/No_Samples&Cor_PathwayActivities_RF.pdf",height = 5, width = 6)
ggplot(df_RF_PA, aes(x = Cor, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) + theme_classic()

dev.off()



#5) Transcription Factors
# Plot number of samples against  results for all the ML methods

#1-1) Ridge
df_Ridge_TF = data.frame(Cor = Ridge_TF[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/Result2/No_Samples&Cor_TranscriptionFactors_Ridge.pdf",height = 5, width = 6)

ggplot(df_Ridge_TF, aes(x = Cor, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE)+ theme_classic()

dev.off()

a = rownames(Ridge_PW_new)[Ridge_PW_new$N_genes!=0]
df_Ridge_TF2 = df_Ridge_TF
rownames(df_Ridge_TF2) = rownames(Ridge_TF)
df_Ridge_TF2 = df_Ridge_TF2[a,]
  
  
cor(df_Ridge_TF2$Cor,df_Ridge_TF2$N)

lm(Cor~N, data=df_Ridge_TF)

#1-2) MLP
df_MLP_TF = data.frame(Cor = MLP_TF[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/Result2/No_Samples&Cor_TranscriptionFactors_MLP.pdf",height = 5, width = 6)

ggplot(df_MLP_TF, aes(x = Cor, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) + theme_classic()

dev.off()

#1-3) Lasso
df_Lasso_TF = data.frame(Cor = Lasso_TF[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/Result2/No_Samples&Cor_TranscriptionFactors_Lasso.pdf",height = 5, width = 6)

ggplot(df_Lasso_TF, aes(x = Cor, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) + theme_classic()

dev.off()

#1-4) ENet
df_ENet_TF = data.frame(Cor = ENet_TF[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/Result2/No_Samples&Cor_TranscriptionFactors_ENet.pdf",height = 5, width = 6)
ggplot(df_ENet_TF, aes(x = Cor, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) + theme_classic()

dev.off()

#1-5) RF
df_RF_TF = data.frame(Cor = RF_TF[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/Result2/No_Samples&Cor_TranscriptionFactors_RF.pdf",height = 5, width = 6)
ggplot(df_RF_TF, aes(x = Cor, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE) + theme_classic()

dev.off()






