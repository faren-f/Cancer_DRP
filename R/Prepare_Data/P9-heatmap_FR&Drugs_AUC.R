rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

sen_PRISM = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")

#Plot Ridge
WholeGenes_Ridge = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF1_WholeGenes_ridge.rds")
Landmark_Ridge = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF2-Landmark_Ridge.rds")
Drug_Pathways_Ridge = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF6-PW_ridge.rds")
TF_activity_Ridge = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF5-decoupleR_gsea2_ridge.rds")
Pathway_activity_Ridge = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF7-progeny_ridge.rds")

Drug_FSMethods_AUC_Ridge = cbind(WholeGenes_Ridge$AUC,
                                     Landmark_Ridge$AUC,
                                     Drug_Pathways_Ridge$AUC,
                                     Pathway_activity_Ridge$AUC,
                                     TF_activity_Ridge$AUC)

rownames(Drug_FSMethods_AUC_Ridge) = colnames(sen_PRISM)
Drug_FSMethods_AUC_Ridge = Drug_FSMethods_AUC_Ridge[-c(3,10,12,14),]  
Drug_FSMethods_AUC_Ridge[Drug_FSMethods_AUC_Ridge<0.5] = 0.5

plt_Ridge = pheatmap::pheatmap(t(Drug_FSMethods_AUC_Ridge),cluster_rows = FALSE, 
                               cluster_cols = FALSE, color = RColorBrewer::brewer.pal(8,"Blues"))

pdf(paste0("Figures/FS/heatmap-FR&Drugs_all_ML/heatmap_FR&Drugs_AUC_Ridge.pdf"), height = 2.8, width = 7)
plt_Ridge
dev.off()


# Plot Lasso
WholeGenes_Lasso = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Lasso/RF1_WholeGenes_Lasso.rds")
Landmark_Lasso = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Lasso/RF2-Landmark_Lasso.rds")
Drug_Pathways_Lasso = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Lasso/RF6-PW_Lasso.rds")
TF_activity_Lasso = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Lasso/RF5-decoupleR_gsea2_Lasso.rds")
Pathway_activity_Lasso = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Lasso/RF7-progeny_Lasso.rds")

Drug_FSMethods_AUC_Lasso = cbind(WholeGenes_Lasso$AUC,
                                     Landmark_Lasso$AUC,
                                     Drug_Pathways_Lasso$AUC,
                                     Pathway_activity_Lasso$AUC,
                                     TF_activity_Lasso$AUC)
rownames(Drug_FSMethods_AUC_Lasso) = colnames(sen_PRISM)
Drug_FSMethods_AUC_Lasso = Drug_FSMethods_AUC_Lasso[-c(3,10,12,14),]  

Drug_FSMethods_AUC_Lasso[Drug_FSMethods_AUC_Lasso<0.5] = 0.5
plt_Lasso = pheatmap::pheatmap(t(Drug_FSMethods_AUC_Lasso),cluster_rows = FALSE, 
                               cluster_cols = FALSE,color = RColorBrewer::brewer.pal(8,"Blues"))
pdf(paste0("Figures/FS/heatmap-FR&Drugs_all_ML/heatmap_FR&Drugs_AUC_Lasso.pdf"), height = 2.8, width = 7)
plt_Lasso
dev.off()


# Plot Elastic Net
WholeGenes_ENet = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ENet/RF1_WholeGenes_ENet.rds")
Landmark_ENet = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ENet/RF2-Landmark_ENet.rds")
Drug_Pathways_ENet = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ENet/RF6-PW_ENet.rds")
TF_activity_ENet = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ENet/RF5-decoupleR_gsea2_ENet.rds")
Pathway_activity_ENet = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ENet/RF7-progeny_ENet.rds")

Drug_FSMethods_AUC_ENet = cbind(WholeGenes_ENet$AUC,
                                    Landmark_ENet$AUC,
                                    Drug_Pathways_ENet$AUC,
                                    Pathway_activity_ENet$AUC,
                                    TF_activity_ENet$AUC)
rownames(Drug_FSMethods_AUC_ENet) = colnames(sen_PRISM)
Drug_FSMethods_AUC_ENet = Drug_FSMethods_AUC_ENet[-c(3,10,12,14),]  

Drug_FSMethods_AUC_ENet[Drug_FSMethods_AUC_ENet<0.5] = 0.5
plt_ENet = pheatmap::pheatmap(t(Drug_FSMethods_AUC_ENet),cluster_rows = FALSE, 
                              cluster_cols = FALSE,color = RColorBrewer::brewer.pal(8,"Blues"))

pdf(paste0("Figures/FS/heatmap-FR&Drugs_all_ML/heatmap_FR&Drugs_AUC_ENet.pdf"), height = 2.8, width = 7)
plt_ENet
dev.off()


# Plot RandomForest
WholeGenes_RF = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/RF/RF1_WholeGenes_RF.rds")
Landmark_RF = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/RF/RF2-Landmark_RF.rds")
Drug_Pathways_RF = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/RF/RF6-PW_RF.rds")
TF_activity_RF = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/RF/RF5-decoupleR_gsea2_RF.rds")
Pathway_activity_RF = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/RF/RF7-progeny_RF.rds")

Drug_FSMethods_AUC_RF = cbind(WholeGenes_RF$AUC,
                                  Landmark_RF$AUC,
                                  Drug_Pathways_RF$AUC,
                                  Pathway_activity_RF$AUC,
                                  TF_activity_RF$AUC)
rownames(Drug_FSMethods_AUC_RF) = colnames(sen_PRISM)
Drug_FSMethods_AUC_RF = Drug_FSMethods_AUC_RF[-c(3,10,12,14),]  

Drug_FSMethods_AUC_RF[Drug_FSMethods_AUC_RF<0.5] = 0.5
plt_RF = pheatmap::pheatmap(t(Drug_FSMethods_AUC_RF),cluster_rows = FALSE, 
                            cluster_cols = FALSE,color = RColorBrewer::brewer.pal(8,"Blues"))

pdf(paste0("Figures/FS/heatmap-FR&Drugs_all_ML/heatmap_FR&Drugs_AUC_RF.pdf"), height = 2.8, width = 7)
plt_RF
dev.off()


# Plot MLP
WholeGenes_MLP = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/MLP/RF1_WholeGenes_MLP.rds")
Landmark_MLP = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/MLP/RF2-Landmark_MLP.rds")
Drug_Pathways_MLP = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/MLP/RF6-PW_MLP.rds")
TF_activity_MLP = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/MLP/RF5-decoupleR_gsea2_MLP.rds")
Pathway_activity_MLP = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/MLP/RF7-progeny_MLP.rds")

Drug_FSMethods_AUC_MLP = cbind(WholeGenes_MLP$AUC,
                                   Landmark_MLP$AUC,
                                   Drug_Pathways_MLP$AUC,
                                   Pathway_activity_MLP$AUC,
                                   TF_activity_MLP$AUC)
rownames(Drug_FSMethods_AUC_MLP) = colnames(sen_PRISM)
Drug_FSMethods_AUC_MLP = Drug_FSMethods_AUC_MLP[-c(3,10,12,14),]  

Drug_FSMethods_AUC_MLP[Drug_FSMethods_AUC_MLP<0.5] = 0.5
plt_MLP = pheatmap::pheatmap(t(Drug_FSMethods_AUC_MLP),cluster_rows = FALSE, 
                             cluster_cols = FALSE, color = RColorBrewer::brewer.pal(8,"Blues"))

pdf(paste0("Figures/FS/heatmap-FR&Drugs_all_ML/heatmap_FR&Drugs_AUC_MLP.pdf"), height = 2.8, width = 7)
plt_MLP
dev.off()


