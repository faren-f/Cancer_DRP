
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
PW = readRDS("Final_Result/Train@PRISM_Test@TCGA_FS/RF6-PW_EachTarget.rds")
PW_ranksum = c()
for(i in names(PW)){
  PW_ranksum = c(PW_ranksum, PW[[i]][[1]][5]) 
}


Drug_FSMethods = cbind(WholeGenes, Landmark, DoRothEA, progeny)
rownames(Drug_FSMethods) = colnames(sen_PRISM)
  
  
Drug_FSMethods_binary = ifelse(Drug_FSMethods<0.05,0,1)
pheatmap::pheatmap(Drug_FSMethods_binary,cluster_rows = FALSE, 
                   cluster_cols = FALSE)
