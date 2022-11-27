
rm(list=ls())
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
sen_PRISM = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")

WholeGenes = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF1_WholeGenes_ridge.rds")
WholeGenes = WholeGenes$Ranksum

Landmark = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF2-Landmark_Ridge.rds")
Landmark = Landmark$Ranksum

DoRothEA = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF3-DoRothEA_ridge.rds")
DoRothEA = DoRothEA$Ranksum

progeny = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF7-progeny_ridge.rds")
progeny = progeny$Ranksum

Drug_Pathways = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF6-PW_ridge.rds")
Drug_Pathways = Drug_Pathways$Ranksum

OncoKB = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF8-OncoKB_ridge.rds")
OncoKB = OncoKB$Ranksum

CancerGenes_GDSC = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF15-CancerGenes_ridge.rds")
CancerGenes_GDSC = CancerGenes_GDSC$Ranksum



# intersect_Landmark_OncoKB = readRDS("Final_Result/Train@PRISM_Test@TCGA_FS/RF13-Landmark_OncoKB_genes.rds")
# intersect_Landmark_OncoKB = intersect_Landmark_OncoKB$Ranksum
# 
# Landmark_OncoKB = readRDS("Final_Result/Train@PRISM_Test@TCGA_FS/RF14-Landmark&OncoKB_genes.rds")
# Landmark_OncoKB = Landmark_OncoKB$Ranksum


# PW_ranksum = c()
# for(i in names(PW)){
#   PW_ranksum = c(PW_ranksum, PW[[i]][[1]][5]) 
# }


Drug_FSMethods = cbind(WholeGenes,Landmark,DoRothEA,progeny,Drug_Pathways,OncoKB,CancerGenes_GDSC)
rownames(Drug_FSMethods) = colnames(sen_PRISM)
Drug_FSMethods = Drug_FSMethods[-c(3,10,12,14),]  
  
Drug_FSMethods_binary = ifelse(Drug_FSMethods<0.05,0,1)
plt = pheatmap::pheatmap(t(Drug_FSMethods_binary),cluster_rows = FALSE, 
                   cluster_cols = FALSE,color = c("deeppink4", "snow2"))


pdf(paste0("Figures/FS/All_Methods/heatmap_FS_Drugs.pdf"), height = 4, width = 5)
plt
dev.off()




