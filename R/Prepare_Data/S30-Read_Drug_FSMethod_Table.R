
rm(list=ls())
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
sen_PRISM = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")

WholeGenes = readRDS("Final_Result/Train@PRISM_Test@TCGA_FS/RF1_WholeGenes.rds")
WholeGenes = WholeGenes$Ranksum
Landmark = readRDS("Final_Result/Train@PRISM_Test@TCGA_FS/RF2-Landmark.rds")
Landmark = Landmark$Ranksum
DoRothEA = readRDS("Final_Result/Train@PRISM_Test@TCGA_FS/RF3-DoRothEA.rds")
DoRothEA = DoRothEA$Ranksum
progeny = readRDS("Final_Result/Train@PRISM_Test@TCGA_FS/RF7-progeny.rds")
progeny = progeny$Ranksum
Drug_Pathways = readRDS("Final_Result/Train@PRISM_Test@TCGA_FS/RF6-PW.rds")
Drug_Pathways = Drug_Pathways$Ranksum

# PW_ranksum = c()
# for(i in names(PW)){
#   PW_ranksum = c(PW_ranksum, PW[[i]][[1]][5]) 
# }

Drug_Pathways = rep(1,24)
Drug_Pathways[c(5,8,11,17,21)] = 0
Drug_FSMethods = cbind(WholeGenes, Landmark, DoRothEA, progeny, Drug_Pathways)
rownames(Drug_FSMethods) = colnames(sen_PRISM)
Drug_FSMethods = Drug_FSMethods[-c(3,10,12,14),]  
  
Drug_FSMethods_binary = ifelse(Drug_FSMethods<0.05,0,1)
plt = pheatmap::pheatmap(Drug_FSMethods_binary,cluster_rows = FALSE, 
                   cluster_cols = FALSE,color = c("deeppink4", "snow2"))

pdf(paste0("Figures/FS/All_Methods/heatmap_FS_Drugs.pdf"), height = 4, width = 5)
plt
dev.off()



