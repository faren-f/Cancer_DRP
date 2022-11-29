rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
sen_PRISM = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")

WholeGenes = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF1_WholeGenes_ridge.rds")
Landmark = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF2-Landmark_Ridge.rds")
Drug_Pathways = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF6-PW_ridge.rds")
TF_activity = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF5-decoupleR_gsea2_ridge.rds")
Pathway_activity = readRDS("Final_Result/TrainPRISM&TestTCGA_FS/Ridge/RF7-progeny_ridge.rds")

#Plot AUC-Ranksum
plot(-log(WholeGenes$Ranksum), WholeGenes$AUC, col = "yellow",
     xlim = c(0,8), ylim = c(0,1), pch = 20, 
     xlab = "-log(p-value)", ylab = "AUC_of_ROC")
par(new = TRUE)
plot(-log(Landmark$Ranksum), Landmark$AUC, col = "red",
     xlim = c(0,8), ylim = c(0,1), pch = 20,
     xlab = "-log(p-value)", ylab = "AUC_of_ROC")
par(new = TRUE)
plot(-log(Drug_Pathways$Ranksum), Drug_Pathways$AUC, col = "green",
     xlim = c(0,8), ylim = c(0,1), pch = 20,
     xlab = "-log(p-value)", ylab = "AUC_of_ROC")
par(new = TRUE)
plot(-log(TF_activity$Ranksum), TF_activity$AUC, col = "blue",
     xlim = c(0,8), ylim = c(0,1), pch = 20,
     xlab = "-log(p-value)", ylab = "AUC_of_ROC")
par(new = TRUE)
plot(-log(Pathway_activity$Ranksum), Pathway_activity$AUC, col = "pink",
     xlim = c(0,8), ylim = c(0,1), pch = 20,
     xlab = "-log(p-value)", ylab = "AUC_of_ROC")

# Box Plot AUC
abline(v = -log(0.05),lty = "dotdash", col = "brown")
boxplot(Result_Landmark$Ranksum, Result_Wholegenes$Ranksum, 
          Result_PW$Ranksum, Result_gsea$Ranksum, Result_progeny$Ranksum,
          names=NA, cex=.1, outline = FALSE)


pdf(paste0("Figures/FS/All_Methods/heatmap_FS_Drugs_"), height = 2.8, width = 7)
plt
dev.off()


