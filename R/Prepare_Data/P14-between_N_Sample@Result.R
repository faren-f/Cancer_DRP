rm(list=ls())
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
sen = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")
sen = sen[,-c(3,10,12,14)]
# Finding Number of samples for each drug
N_Na = data.frame(apply(sen,2,function(x){sum(is.na(x))}))
N_S = data.frame(475-N_Na[,1])
rownames(N_S) = rownames(N_Na)


# No samples/AUC, No samples/-log(Ranksum) for Whole Genes  
WholeGenes = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF1_WholeGenes_ridge.rds")
WholeGenes = WholeGenes$AUC
#WholeGenes = -log(WholeGenes)
WholeGenes = WholeGenes[-c(3,10,12,14)]

df_Ridge_WG = data.frame(Cor = WholeGenes, N = N_S[,1])
ggplot(df_Ridge_WG, aes(x = Cor, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE)+ theme_classic()


Landmark = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF2-Landmark_Ridge.rds")
Landmark = Landmark$Ranksum
Landmark = Landmark[-c(3,10,12,14)]

df_Ridge_Landmark = data.frame(Ranksum = Landmark, N = N_S[,1])
ggplot(df_Ridge_Landmark, aes(x = Ranksum, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE)+ theme_classic()


TF = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF5-decoupleR_gsea2_ridge.rds")
TF = TF$Ranksum
TF = TF[-c(3,10,12,14)]

df_Ridge_TF = data.frame(wilkoxtest = TF, N = N_S[,1])
ggplot(df_Ridge_TF, aes(x = wilkoxtest, y = N)) + geom_point(size = 0.5) +
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE)+ theme_classic()


progeny = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF7-progeny_ridge.rds")
progeny = progeny$Ranksum

Drug_Pathways = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF6-PW_ridge.rds")
Drug_Pathways = Drug_Pathways$Ranksum














