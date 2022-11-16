rm(list=ls())
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

sen_PRISM = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")

WholeGenes = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ridge/RF1_WholeGenes_ridge.rds")
WholeGenes = WholeGenes$Ranksum

Landmark = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF2-Landmark_Ridge.rds")
Landmark = Landmark$Ranksum

Drug_Pathways = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ridge/RF6-PW2_ridge.rds")
Drug_Pathways = Drug_Pathways$Ranksum

#OncoKB = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ENet/)
#OncoKB = OncoKB$Ranksum

# DoRothEA = readRDS("Final_Result/Train@PRISM_Test@TCGA_FS/RF3-DoRothEA.rds")
# DoRothEA = DoRothEA$Ranksum

progeny = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ridge/RF7-progeny_ridge.rds")
progeny = progeny$Ranksum

dR_gsea = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ridge/RF5-decoupleR_gsea2_ridge.rds")
dR_gsea = dR_gsea$Ranksum

# CancerGenes_GDSC = readRDS("Final_Result/Train@PRISM_Test@TCGA_FS/RF15-CancerGenes.rds")
# CancerGenes_GDSC = CancerGenes_GDSC$Ranksum
# 
# 
# dR_aucell = readRDS("Final_Result/Train@PRISM_Test@TCGA_FS/RF5-decoupleR_aucell.rds")
# dR_aucell = dR_aucell$Ranksum
# 
# dR_consensus = readRDS("Final_Result/Train@PRISM_Test@TCGA_FS/RF5-decoupleR_consensus.rds")
# dR_consensus = dR_consensus$Ranksum
# 
# dR_gsva = readRDS("Final_Result/Train@PRISM_Test@TCGA_FS/RF5-decoupleR_gsva.rds")
# dR_gsva = dR_gsva$Ranksum
# 
# dR_mdt = readRDS("Final_Result/Train@PRISM_Test@TCGA_FS/RF5-decoupleR_mdt.rds")
# dR_mdt = dR_mdt$Ranksum
# 
# dR_mlm = readRDS("Final_Result/Train@PRISM_Test@TCGA_FS/RF5-decoupleR_mlm.rds")
# dR_mlm = dR_mlm$Ranksum
# 
# dR_udt = readRDS("Final_Result/Train@PRISM_Test@TCGA_FS/RF5-decoupleR_udt.rds")
# dR_udt = dR_udt$Ranksum
# 
# dR_ulm = readRDS("Final_Result/Train@PRISM_Test@TCGA_FS/RF5-decoupleR_ulm.rds")
# dR_ulm = dR_ulm$Ranksum
# 
# dR_viper = readRDS("Final_Result/Train@PRISM_Test@TCGA_FS/RF5-decoupleR_viper.rds")
# dR_viper = dR_viper$Ranksum
# 
# dR_gsea = readRDS("Final_Result/Train@PRISM_Test@TCGA_FS/RF5-decoupleR_gsea2.rds")
# dR_gsea = dR_gsea$Ranksum

#dR_ora = readRDS("Final_Result/Train@PRISM_Test@TCGA_FS/RF5")
#dR_ora = dR_ora$Ranksum

# dR_wmean = readRDS("Final_Result/Train@PRISM_Test@TCGA_FS/RF5-decoupleR_wmean2.rds")
# dR_wmean = dR_wmean$Ranksum
# 
# dR_wsum = readRDS("Final_Result/Train@PRISM_Test@TCGA_FS/RF5-decoupleR_wsum2.rds")
# dR_wsum = dR_wsum$Ranksum


# Drug_FSMethods = cbind(WholeGenes,Landmark,DoRothEA,progeny,Drug_Pathways,OncoKB,CancerGenes_GDSC,
#                        dR_aucell,dR_consensus,dR_gsva,dR_mdt,dR_mlm,dR_udt,dR_ulm,
#                        dR_viper, dR_gsea,dR_wmean,dR_wsum)
TF_activity = dR_gsea
Pathway_activity = progeny
Drug_FSMethods = cbind(WholeGenes, Landmark, Drug_Pathways, TF_activity, Pathway_activity)

rownames(Drug_FSMethods) = colnames(sen_PRISM)
Drug_FSMethods = Drug_FSMethods[-c(3,10,12,14),]  

Drug_FSMethods_binary = ifelse(Drug_FSMethods<0.05,0,1)
plt = pheatmap::pheatmap(t(Drug_FSMethods_binary),cluster_rows = FALSE, 
                         cluster_cols = FALSE,color = c("deeppink4", "snow2"))


pdf(paste0("Figures/FS/All_Methods/heatmap_FR_ridge.pdf"), height = 4, width = 5)
plt
dev.off()




