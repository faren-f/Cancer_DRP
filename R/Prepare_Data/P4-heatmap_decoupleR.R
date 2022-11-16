rm(list=ls())
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
sen_PRISM = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")


dR_aucell = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ridge/RF5-decoupleR_aucell_ridge.rds")
dR_aucell = dR_aucell$Ranksum

dR_consensus = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ridge/RF5-decoupleR_consensus_ridge.rds")
dR_consensus = dR_consensus$Ranksum

dR_gsva = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ridge/RF5-decoupleR_gsva_ridge.rds")
dR_gsva = dR_gsva$Ranksum

dR_mdt = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ridge/RF5-decoupleR_mdt_ridge.rds")
dR_mdt = dR_mdt$Ranksum

dR_mlm = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ridge/RF5-decoupleR_mlm_ridge.rds")
dR_mlm = dR_mlm$Ranksum

dR_udt = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ridge/RF5-decoupleR_udt_ridge.rds")
dR_udt = dR_udt$Ranksum

dR_ulm = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ridge/RF5-decoupleR_ulm_ridge.rds")
dR_ulm = dR_ulm$Ranksum

dR_viper = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ridge/RF5-decoupleR_viper_ridge.rds")
dR_viper = dR_viper$Ranksum

dR_gsea = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ridge/RF5-decoupleR_gsea2_ridge.rds")
dR_gsea = dR_gsea$Ranksum

#dR_ora = readRDS("Final_Result/Train@PRISM_Test@TCGA_FS/RF5")
#dR_ora = dR_ora$Ranksum

dR_wmean = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ridge/RF5-decoupleR_wmean2_ridge.rds")
dR_wmean = dR_wmean$Ranksum

dR_wsum = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/ridge/RF5-decoupleR_wsum2_ridge.rds")
dR_wsum = dR_wsum$Ranksum



Drug_FSMethods = cbind(dR_aucell,dR_consensus,dR_gsva,dR_mdt,dR_mlm,dR_udt,dR_ulm,
                       dR_viper,dR_gsea,dR_wmean,dR_wsum)
rownames(Drug_FSMethods) = colnames(sen_PRISM)
Drug_FSMethods = Drug_FSMethods[-c(3,10,12,14),]  

Drug_FSMethods_binary = ifelse(Drug_FSMethods<0.05,0,1)
plt = pheatmap::pheatmap(t(Drug_FSMethods_binary),cluster_rows = FALSE, 
                         cluster_cols = FALSE,color = c("deeppink4", "snow2"))


pdf(paste0("Figures/FS/DecoupleR/heatmap_decoupleR_ridge.pdf"), height = 4, width = 5)
plt
dev.off()




