rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
sen_PRISM = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")

WholeGenes_Ridge = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF1_WholeGenes_ridge.rds")[-c(3,10,12,14),]
Landmark_Ridge = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF2-Landmark_Ridge.rds")[-c(3,10,12,14),]
Drug_Pathways_Ridge = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF6-PW_ridge.rds")[-c(3,10,12,14),]
TF_activity_Ridge = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF5-decoupleR_gsea2_ridge.rds")[-c(3,10,12,14),]
Pathway_activity_Ridge = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF7-progeny_ridge.rds")[-c(3,10,12,14),]

WholeGenes_Lasso = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Lasso/RF1_WholeGenes_Lasso.rds")[-c(3,10,12,14),]
Landmark_Lasso = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Lasso/RF2-Landmark_Lasso.rds")[-c(3,10,12,14),]
Drug_Pathways_Lasso = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Lasso/RF6-PW_Lasso.rds")[-c(3,10,12,14),]
TF_activity_Lasso = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Lasso/RF5-decoupleR_gsea2_Lasso.rds")[-c(3,10,12,14),]
Pathway_activity_Lasso = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Lasso/RF7-progeny_Lasso.rds")[-c(3,10,12,14),]

WholeGenes_ENet = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ENet/RF1_WholeGenes_ENet.rds")[-c(3,10,12,14),]
Landmark_ENet = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ENet/RF2-Landmark_ENet.rds")[-c(3,10,12,14),]
Drug_Pathways_ENet = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ENet/RF6-PW_ENet.rds")[-c(3,10,12,14),]
TF_activity_ENet = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ENet/RF5-decoupleR_gsea2_ENet.rds")[-c(3,10,12,14),]
Pathway_activity_ENet = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ENet/RF7-progeny_ENet.rds")[-c(3,10,12,14),]

WholeGenes_RF = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/RF/RF1_WholeGenes_RF.rds")[-c(3,10,12,14),]
Landmark_RF = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/RF/RF2-Landmark_RF.rds")[-c(3,10,12,14),]
Drug_Pathways_RF = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/RF/RF6-PW_RF.rds")[-c(3,10,12,14),]
TF_activity_RF = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/RF/RF5-decoupleR_gsea2_RF.rds")[-c(3,10,12,14),]
Pathway_activity_RF = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/RF/RF7-progeny_RF.rds")[-c(3,10,12,14),]

WholeGenes_MLP = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/MLP/RF1_WholeGenes_MLP.rds")[-c(3,10,12,14),]
Landmark_MLP = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/MLP/RF2-Landmark_MLP.rds")[-c(3,10,12,14),]
Drug_Pathways_MLP = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/MLP/RF6-PW_MLP.rds")[-c(3,10,12,14),]
TF_activity_MLP = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/MLP/RF5-decoupleR_gsea2_MLP.rds")[-c(3,10,12,14),]
Pathway_activity_MLP = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/MLP/RF7-progeny_MLP.rds")[-c(3,10,12,14),]


#Grouped barplot for Ranksum Test
WholeGenes_Ridge_ranksum = sum(WholeGenes_Ridge$Ranksum<0.05)
Landmark_Ridge_ranksum = sum(Landmark_Ridge$Ranksum<0.05)
Drug_Pathways_Ridge_ranksum = sum(Drug_Pathways_Ridge$Ranksum<0.05)
TF_activity_Ridge_ranksum = sum(TF_activity_Ridge$Ranksum<0.05)
Pathway_activity_Ridge_ranksum = sum(Pathway_activity_Ridge$Ranksum<0.05)

WholeGenes_Lasso_ranksum = sum(WholeGenes_Lasso$Ranksum<0.05)
Landmark_Lasso_ranksum = sum(Landmark_Lasso$Ranksum<0.05)
Drug_Pathways_Lasso_ranksum = sum(Drug_Pathways_Lasso$Ranksum<0.05)
TF_activity_Lasso_ranksum = sum(TF_activity_Lasso$Ranksum<0.05)
Pathway_activity_Lasso_ranksum = sum(Pathway_activity_Lasso$Ranksum<0.05)

WholeGenes_ENet_ranksum = sum(WholeGenes_ENet$Ranksum<0.05)
Landmark_ENet_ranksum = sum(Landmark_ENet$Ranksum<0.05)
Drug_Pathways_ENet_ranksum = sum(Drug_Pathways_ENet$Ranksum<0.05)
TF_activity_ENet_ranksum = sum(TF_activity_ENet$Ranksum<0.05)
Pathway_activity_ENet_ranksum = sum(Pathway_activity_ENet$Ranksum<0.05)

WholeGenes_RF_ranksum = sum(WholeGenes_RF$Ranksum<0.05)
Landmark_RF_ranksum = sum(Landmark_RF$Ranksum<0.05)
Drug_Pathways_RF_ranksum = sum(Drug_Pathways_RF$Ranksum<0.05)
TF_activity_RF_ranksum = sum(TF_activity_RF$Ranksum<0.05)
Pathway_activity_RF_ranksum = sum(Pathway_activity_RF$Ranksum<0.05)

WholeGenes_MLP_ranksum = sum(WholeGenes_MLP$Ranksum<0.05)
Landmark_MLP_ranksum = sum(Landmark_MLP$Ranksum<0.05)
Drug_Pathways_MLP_ranksum = sum(Drug_Pathways_MLP$Ranksum<0.05)
TF_activity_MLP_ranksum = sum(TF_activity_MLP$Ranksum<0.05)
Pathway_activity_MLP_ranksum = sum(Pathway_activity_MLP$Ranksum<0.05)

ML_methods = factor(rep(c("Ridge","Lasso","ENet","RF","MLP"),5),
                       levels = c("Ridge","MLP","Lasso","ENet","RF"))
FR_methods = factor(rep(c("WG","LM","PW","PA","TF"),c(5,5,5,5,5)),
                        levels = c("WG","LM","PW","PA","TF"))
N_sig_drugs = c(WholeGenes_Ridge_ranksum,WholeGenes_Lasso_ranksum,
  WholeGenes_ENet_ranksum,WholeGenes_RF_ranksum,WholeGenes_MLP_ranksum,
  Landmark_Ridge_ranksum,Landmark_Lasso_ranksum,Landmark_ENet_ranksum,
  Landmark_RF_ranksum,Landmark_MLP_ranksum,Drug_Pathways_Ridge_ranksum,
  Drug_Pathways_Lasso_ranksum,Drug_Pathways_ENet_ranksum,Drug_Pathways_RF_ranksum,
  Drug_Pathways_MLP_ranksum,Pathway_activity_Ridge_ranksum,Pathway_activity_Lasso_ranksum,
  Pathway_activity_ENet_ranksum,Pathway_activity_RF_ranksum,Pathway_activity_MLP_ranksum,
  TF_activity_Ridge_ranksum,TF_activity_Lasso_ranksum,TF_activity_ENet_ranksum,
  TF_activity_RF_ranksum,TF_activity_MLP_ranksum)

data = data.frame(FR_methods,ML_methods,N_sig_drugs)
ggplot(data, aes(fill=FR_methods, y=N_sig_drugs, x=ML_methods)) + 
  geom_bar(position="dodge", stat="identity",width = 0.5, color = "black", size = 0.2)+
  scale_fill_manual(values=c("#F2E6D6","#A0E0E4","#FCC0C5","#A49393","#F582A8")) + 
  theme_minimal()


#Grouped barplot for AUC
WholeGenes_Ridge_AUC = mean(WholeGenes_Ridge$AUC)
Landmark_Ridge_AUC = mean(Landmark_Ridge$AUC)
Drug_Pathways_Ridge_AUC = mean(Drug_Pathways_Ridge$AUC)
TF_activity_Ridge_AUC = mean(TF_activity_Ridge$AUC)
Pathway_activity_Ridge_AUC = mean(Pathway_activity_Ridge$AUC)

WholeGenes_Lasso_AUC = mean(WholeGenes_Lasso$AUC)
Landmark_Lasso_AUC = mean(Landmark_Lasso$AUC)
Drug_Pathways_Lasso_AUC = mean(Drug_Pathways_Lasso$AUC)
TF_activity_Lasso_AUC = mean(TF_activity_Lasso$AUC)
Pathway_activity_Lasso_AUC = mean(Pathway_activity_Lasso$AUC)

WholeGenes_ENet_AUC = mean(WholeGenes_ENet$AUC)
Landmark_ENet_AUC = mean(Landmark_ENet$AUC)
Drug_Pathways_ENet_AUC = mean(Drug_Pathways_ENet$AUC)
TF_activity_ENet_AUC = mean(TF_activity_ENet$AUC)
Pathway_activity_ENet_AUC = mean(Pathway_activity_ENet$AUC)

WholeGenes_RF_AUC = mean(WholeGenes_RF$AUC)
Landmark_RF_AUC = mean(Landmark_RF$AUC)
Drug_Pathways_RF_AUC = mean(Drug_Pathways_RF$AUC)
TF_activity_RF_AUC = mean(TF_activity_RF$AUC)
Pathway_activity_RF_AUC = mean(Pathway_activity_RF$AUC)

WholeGenes_MLP_AUC = mean(WholeGenes_MLP$AUC)
Landmark_MLP_AUC = mean(Landmark_MLP$AUC)
Drug_Pathways_MLP_AUC = mean(Drug_Pathways_MLP$AUC)
TF_activity_MLP_AUC = mean(TF_activity_MLP$AUC)
Pathway_activity_MLP_AUC = mean(Pathway_activity_MLP$AUC)

ML_methods = factor(rep(c("Ridge","Lasso","ENet","RF","MLP"),5),
                    levels = c("Ridge","MLP","Lasso","ENet","RF"))
FR_methods = factor(rep(c("WG","LM","PW","PA","TF"),c(5,5,5,5,5)),
                    levels = c("WG","LM","PW","PA","TF"))
N_sig_drugs = c(WholeGenes_Ridge_AUC,WholeGenes_Lasso_AUC,
                WholeGenes_ENet_AUC,WholeGenes_RF_AUC,WholeGenes_MLP_AUC,
                Landmark_Ridge_AUC,Landmark_Lasso_AUC,Landmark_ENet_AUC,
                Landmark_RF_AUC,Landmark_MLP_AUC,Drug_Pathways_Ridge_AUC,
                Drug_Pathways_Lasso_AUC,Drug_Pathways_ENet_AUC,Drug_Pathways_RF_AUC,
                Drug_Pathways_MLP_AUC,Pathway_activity_Ridge_AUC,Pathway_activity_Lasso_AUC,
                Pathway_activity_ENet_AUC,Pathway_activity_RF_AUC,Pathway_activity_MLP_AUC,
                TF_activity_Ridge_AUC,TF_activity_Lasso_AUC,TF_activity_ENet_AUC,
                TF_activity_RF_AUC,TF_activity_MLP_AUC)

data = data.frame(FR_methods,ML_methods,N_sig_drugs)
ggplot(data, aes(fill=FR_methods, y=N_sig_drugs, x=ML_methods)) + 
  geom_bar(position="dodge", stat="identity",width = 0.5, color = "black", size = 0.2)+
  scale_fill_manual(values=c("#F2E6D6","#A0E0E4","#FCC0C5","#A49393","#F582A8")) + 
  theme_minimal()

