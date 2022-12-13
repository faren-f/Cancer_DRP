rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
sen = readRDS("Processed_data/S1/sensitivity_matrix_AUC.rds")
GE = readRDS("Processed_data/S1/expresion_matrix.rds")
N_drugs = 1448

#Landmark
l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")

#TF activity
TF = read.table("Processed_from_Python/TF_gsea2_PRISM/TF_gsea2_PRISM.csv",
                sep = ",",header = TRUE, row.names = 1)

# Progeny
PA = readRDS("Processed_data/Other/pw_act_GE.rds")

#PW
#Read Drug Pathway results
Ridge_PW = c()
I_zeros = c()
for(i in 1:N_drugs){
  print(i)
  R_PW = readRDS(paste0("Processed_from_SLURM/Results_DrugPathways_All_Models/Result_",as.character(i),".rds"))
  if(is.null(nrow(R_PW))){
    I_zeros = c(I_zeros,i)
  }else{
    Ridge_PW = rbind(Ridge_PW, R_PW[4,])
  }
}


pdf(paste0("Figures/FS/Number_of_genes_in_FR_Methods/BoxPlot_N_genes_FR_Methods.pdf"),
    height = 5, width = 7)
boxplot(ncol(GE), length(l1000_genes), Ridge_PW[,5],
        ncol(PA),
        ncol(TF),
        boxfill = "#FCC0C5")  # Box fill color)
#col = c("#F2E6D6","#A0E0E4","#FCC0C5","#A49393","#F582A8"),
#medcol=c("#F2E6D6","#A0E0E4","#FCC0C5","#A49393","#F582A8"),
#border = c("#F2E6D6","#A0E0E4","#FCC0C5","#A49393","#F582A8"))
# boxcol = 2   # Box border color

dev.off()


pdf(paste0("Figures/FS/Number_of_genes_in_FR_Methods/BarPlot_N_genes_FR_Methods.pdf"),
    height = 5, width = 7)

data = data.frame(
  name=letters[1:5],
  value = c(ncol(GE), length(l1000_genes), mean(Ridge_PW[,5]),
            ncol(PA), ncol(TF)),
  sd=c(0, 0, sd(Ridge_PW[,5]), 0, 0))

ggplot(data) +
  geom_bar(aes(x=name, y=value), stat="identity", 
           fill= c("#F2E6D6","#A0E0E4","#FCC0C5","#A49393","#F582A8"), alpha=1, 
           width = 0.7) +
  geom_errorbar(aes(x=name, ymin=value-sd, ymax=value+sd), width=0.3, 
                colour="black", alpha = 0.8, size = 0.9)+ theme_classic()

dev.off()

