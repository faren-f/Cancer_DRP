rm(list=ls())

library(ggplot2)
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
sen = readRDS("Processed_data/S1/sensitivity_matrix_AUC.rds")
N_Na = data.frame(apply(sen,2,function(x){sum(is.na(x))}))
N_S = data.frame(475-N_Na[,1])

N_drugs = 1448
RF_Landmark = c()
ENet_Landmark = c()
Lasso_Landmark = c()
Ridge_Landmark = c()
MLP_Landmark = c()

for(i in 1:N_drugs){
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


# hist(RF[,1],30)
# hist(Ridge[,1],30)
# hist(ENet[,1],30)
# hist(Lasso[,1],30)
# hist(MLP[,1],30)


pdf("Figures/FS/Results_whithin/MeanCorr_Landmark_All_ML.pdf", height = 5, width = 7)
boxplot(cbind(Ridge_Landmark[,1], MLP_Landmark[,1], 
              Lasso_Landmark[,1], ENet_Landmark[,1], 
              RF_Landmark[,1]), names=NA, cex=.5, 
        main="MeanCorr", col = "#A0E0E4")

dev.off()

pdf("Figures/FS/Results_whithin/STDCorr_Landmark_All_ML.pdf", height = 5, width = 7)
boxplot(cbind(Ridge_Landmark[,2], MLP_Landmark[,2], 
              Lasso_Landmark[,2], ENet_Landmark[,2], 
              RF_Landmark[,2]), names=NA, cex=.5, 
        main="STDCorr", col = "#A0E0E4")

dev.off()

pdf("Figures/FS/Results_whithin/MeanMSE_Landmark_All_ML.pdf", height = 5, width = 7)
boxplot(cbind(Ridge_Landmark[,3], MLP_Landmark[,3], 
              Lasso_Landmark[,3], ENet_Landmark[,3], 
              RF_Landmark[,3]), names=NA, cex=.5, 
        main="MeanMSE", col = "#A0E0E4")

dev.off()

pdf("Figures/FS/Results_whithin/STDMSE_Landmark_All_ML.pdf", height = 5, width = 7)
boxplot(cbind(Ridge_Landmark[,4], MLP_Landmark[,4], 
              Lasso_Landmark[,4], ENet_Landmark[,4], 
              RF_Landmark[,4]), names=NA, cex=.5, 
        main="STDMSE", col = "#A0E0E4")

dev.off()

mean(RF[,1])
mean(Ridge[,1])
mean(ENet[,1])
mean(Lasso[,1])
mean(MLP[,1])


# How number of samples related to results?
df_Ridge = data.frame(Cor = Ridge_Landmark[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/No_Samples&Cor_Landmark_Ridge.pdf",height = 5, width = 6)
ggplot(df_Ridge, aes(x = Cor, y = N)) + 
  geom_point(size = 0.5) + geom_smooth(formula = y ~ x, method = "lm", se = FALSE)+
  theme_classic()
dev.off()


df_MLP = data.frame(Cor = MLP_Landmark[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/No_Samples&Cor_Landmark_MLP.pdf",height = 5, width = 6)
ggplot(df_MLP, aes(x = Cor, y = N)) + 
  geom_point(size = 0.5) + geom_smooth(formula = y ~ x, method = "lm", se = FALSE) +
  theme_classic()
dev.off()


df_Lasso = data.frame(Cor = Lasso_Landmark[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/No_Samples&Cor_Landmark_Lasso.pdf",height = 5, width = 6)
ggplot(df_Lasso, aes(x = Cor, y = N)) + 
  geom_point(size = 0.5) + geom_smooth(formula = y ~ x, method = "lm", se = FALSE) +
  theme_classic()
dev.off()


df_ENet = data.frame(Cor = ENet_Landmark[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/No_Samples&Cor_Landmark_ENet.pdf",height = 5, width = 6)
ggplot(df_ENet, aes(x = Cor, y = N)) + 
  geom_point(size = 0.5) + geom_smooth(formula = y ~ x, method = "lm", se = FALSE) +
  theme_classic()
dev.off()


df_RF = data.frame(Cor = RF_Landmark[,1], N = N_S[,1])

pdf("Figures/FS/Results_whithin/No_Samples&Cor_Landmark_RF.pdf",height = 5, width = 6)
ggplot(df_RF, aes(x = Cor, y = N)) + 
  geom_point(size = 0.5) + geom_smooth(formula = y ~ x, method = "lm", se = FALSE) +
  theme_classic()
dev.off()

