rm(list=ls())
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

GE = readRDS("Processed_Data/Step1/expresion_matrix.rds")
sen = readRDS("Processed_Data/Step1/sensitivity_matrix.rds")
dim(GE)
dim(sen)
Rep = 1

mean_mse = rep(0,ncol(sen))
sd_mse = rep(0,ncol(sen))
mean_corr = rep(0,ncol(sen))
sd_corr = rep(0,ncol(sen))
for (i in ncol(sen)){
  
  GE_i = GE[!is.na(sen[,i]),]
  dim(GE_i)
  sen_i = sen_i[!is.na(sen[,i])]
  length(sen_i)
  
  mse = rep(0,Rep)
  corr = rep(0,Rep)
  for(j in 1:Rep){
    
    
    
    
    
  }
  
  
}