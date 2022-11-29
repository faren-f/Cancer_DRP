rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

sen_PRISM = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")

#Plot Ridge
WholeGenes_Ridge = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF1_WholeGenes_ridge.rds")
Landmark_Ridge = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF2-Landmark_Ridge.rds")
Drug_Pathways_Ridge = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF6-PW_ridge.rds")
TF_activity_Ridge = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF5-decoupleR_gsea2_ridge.rds")
Pathway_activity_Ridge = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF7-progeny_ridge.rds")


Drug_FSMethods_ranksum_Ridge = cbind(WholeGenes_Ridge$Ranksum,
                               Landmark_Ridge$Ranksum,
                               Drug_Pathways_Ridge$Ranksum,
                               Pathway_activity_Ridge$Ranksum,
                               TF_activity_Ridge$Ranksum)
rownames(Drug_FSMethods_ranksum_Ridge) = colnames(sen_PRISM)
Drug_FSMethods_ranksum_Ridge = Drug_FSMethods_ranksum_Ridge[-c(3,10,12,14),]  

Drug_FSMethods_ranksum_Ridge = ifelse(Drug_FSMethods_ranksum_Ridge<0.05,0,1)
plt_Ridge = pheatmap::pheatmap(t(Drug_FSMethods_ranksum_Ridge),cluster_rows = FALSE, 
                         cluster_cols = FALSE,color = c("red","white"))

pdf(paste0("Figures/FS/heatmap-FR&Drugs_all_ML/heatmap_FR&Drugs_Ranksum_Ridge.pdf"), height = 2.8, width = 7)
plt_Ridge
dev.off()


# Plot Lasso
WholeGenes_Lasso = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Lasso/RF1_WholeGenes_Lasso.rds")
Landmark_Lasso = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Lasso/RF2-Landmark_Lasso.rds")
Drug_Pathways_Lasso = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Lasso/RF6-PW_Lasso.rds")
TF_activity_Lasso = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Lasso/RF5-decoupleR_gsea2_Lasso.rds")
Pathway_activity_Lasso = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Lasso/RF7-progeny_Lasso.rds")

Drug_FSMethods_ranksum_Lasso = cbind(WholeGenes_Lasso$Ranksum,
                                     Landmark_Lasso$Ranksum,
                                     Drug_Pathways_Lasso$Ranksum,
                                     Pathway_activity_Lasso$Ranksum,
                                     TF_activity_Lasso$Ranksum)
rownames(Drug_FSMethods_ranksum_Lasso) = colnames(sen_PRISM)
Drug_FSMethods_ranksum_Lasso = Drug_FSMethods_ranksum_Lasso[-c(3,10,12,14),]  

Drug_FSMethods_ranksum_Lasso = ifelse(Drug_FSMethods_ranksum_Lasso<0.05,0,1)
plt_Lasso = pheatmap::pheatmap(t(Drug_FSMethods_ranksum_Lasso),cluster_rows = FALSE, 
                               cluster_cols = FALSE,color = c("red","white"))
pdf(paste0("Figures/FS/heatmap-FR&Drugs_all_ML/heatmap_FR&Drugs_Ranksum_Lasso.pdf"), height = 2.8, width = 7)
plt_Lasso
dev.off()


# Plot Elastic Net
WholeGenes_ENet = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ENet/RF1_WholeGenes_ENet.rds")
Landmark_ENet = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ENet/RF2-Landmark_ENet.rds")
Drug_Pathways_ENet = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ENet/RF6-PW_ENet.rds")
TF_activity_ENet = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ENet/RF5-decoupleR_gsea2_ENet.rds")
Pathway_activity_ENet = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ENet/RF7-progeny_ENet.rds")

Drug_FSMethods_ranksum_ENet = cbind(WholeGenes_ENet$Ranksum,
                                     Landmark_ENet$Ranksum,
                                     Drug_Pathways_ENet$Ranksum,
                                     Pathway_activity_ENet$Ranksum,
                                     TF_activity_ENet$Ranksum)
rownames(Drug_FSMethods_ranksum_ENet) = colnames(sen_PRISM)
Drug_FSMethods_ranksum_ENet = Drug_FSMethods_ranksum_ENet[-c(3,10,12,14),]  

Drug_FSMethods_ranksum_ENet = ifelse(Drug_FSMethods_ranksum_ENet<0.05,0,1)
plt_ENet = pheatmap::pheatmap(t(Drug_FSMethods_ranksum_ENet),cluster_rows = FALSE, 
                               cluster_cols = FALSE,color = c("red","white"))

pdf(paste0("Figures/FS/heatmap-FR&Drugs_all_ML/heatmap_FR&Drugs_Ranksum_ENet.pdf"), height = 2.8, width = 7)
plt_ENet
dev.off()


# Plot RandomForest
WholeGenes_RF = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/RF/RF1_WholeGenes_RF.rds")
Landmark_RF = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/RF/RF2-Landmark_RF.rds")
Drug_Pathways_RF = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/RF/RF6-PW_RF.rds")
TF_activity_RF = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/RF/RF5-decoupleR_gsea2_RF.rds")
Pathway_activity_RF = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/RF/RF7-progeny_RF.rds")

Drug_FSMethods_ranksum_RF = cbind(WholeGenes_RF$Ranksum,
                                    Landmark_RF$Ranksum,
                                    Drug_Pathways_RF$Ranksum,
                                    Pathway_activity_RF$Ranksum,
                                    TF_activity_RF$Ranksum)
rownames(Drug_FSMethods_ranksum_RF) = colnames(sen_PRISM)
Drug_FSMethods_ranksum_RF = Drug_FSMethods_ranksum_RF[-c(3,10,12,14),]  

Drug_FSMethods_ranksum_RF = ifelse(Drug_FSMethods_ranksum_RF<0.05,0,1)
plt_RF = pheatmap::pheatmap(t(Drug_FSMethods_ranksum_RF),cluster_rows = FALSE, 
                              cluster_cols = FALSE,color = c("red","white"))

pdf(paste0("Figures/FS/heatmap-FR&Drugs_all_ML/heatmap_FR&Drugs_Ranksum_RF.pdf"), height = 2.8, width = 7)
plt_RF
dev.off()


# Plot MLP
WholeGenes_MLP = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/MLP/RF1_WholeGenes_MLP.rds")
Landmark_MLP = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/MLP/RF2-Landmark_MLP.rds")
Drug_Pathways_MLP = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/MLP/RF6-PW_MLP.rds")
TF_activity_MLP = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/MLP/RF5-decoupleR_gsea2_MLP.rds")
Pathway_activity_MLP = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/MLP/RF7-progeny_MLP.rds")

Drug_FSMethods_ranksum_MLP = cbind(WholeGenes_MLP$Ranksum,
                                  Landmark_MLP$Ranksum,
                                  Drug_Pathways_MLP$Ranksum,
                                  Pathway_activity_MLP$Ranksum,
                                  TF_activity_MLP$Ranksum)
rownames(Drug_FSMethods_ranksum_MLP) = colnames(sen_PRISM)
Drug_FSMethods_ranksum_MLP = Drug_FSMethods_ranksum_MLP[-c(3,10,12,14),]  

Drug_FSMethods_ranksum_MLP = ifelse(Drug_FSMethods_ranksum_MLP<0.05,0,1)
plt_MLP = pheatmap::pheatmap(t(Drug_FSMethods_ranksum_MLP),cluster_rows = FALSE, 
                            cluster_cols = FALSE,color = c("red","white"))

pdf(paste0("Figures/FS/heatmap-FR&Drugs_all_ML/heatmap_FR&Drugs_Ranksum_MLP.pdf"), height = 2.8, width = 7)
plt_MLP
dev.off()



