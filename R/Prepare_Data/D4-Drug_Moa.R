rm(list = ls())
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
sen = readRDS("Processed_data/S1/sensitivity_matrix_AUC.rds")
response = read.csv("Raw_data/PRISM/Secondary/secondary-screen-dose-response-curve-parameters.csv")

N_drugs = 1448
RF_Landmark = c()
ENet_Landmark = c()
Lasso_Landmark = c()
Ridge_Landmark = c()
MLP_Landmark = c()

for(i in 1:N_drugs){
  print(i)
  R_Landmark = readRDS(paste0("Processed_from_SLURM/Results_Landmark_All_Models/Result_",as.character(i),".rds"))
  
  Ridge_Landmark = rbind(Ridge_Landmark, R_Landmark[4,])
  MLP_Landmark = rbind(MLP_Landmark, R_Landmark[5,])
  Lasso_Landmark = rbind(Lasso_Landmark, R_Landmark[3,])
  ENet_Landmark = rbind(ENet_Landmark, R_Landmark[2,])
  RF_Landmark = rbind(RF_Landmark, R_Landmark[1,])
}

rownames(Ridge_Landmark) = colnames(sen)
rownames(MLP_Landmark) = colnames(sen)
rownames(Lasso_Landmark) = colnames(sen)
rownames(ENet_Landmark) = colnames(sen)
rownames(RF_Landmark) = colnames(sen)


#Read Drug Pathway results
RF_PW = c()
ENet_PW = c()
Lasso_PW = c()
Ridge_PW = c()
MLP_PW = c()
I_zeros = c()
for(i in 1:N_drugs){
  print(i)
  R_PW = readRDS(paste0("Processed_from_SLURM/Results_DrugPathways_All_Models/Result_",as.character(i),".rds"))
  if(is.null(nrow(R_PW))){
    I_zeros = c(I_zeros,i)
  }else{
    Ridge_PW = rbind(Ridge_PW, R_PW[4,])
    MLP_PW = rbind(MLP_PW, R_PW[5,])
    Lasso_PW = rbind(Lasso_PW, R_PW[3,])
    ENet_PW = rbind(ENet_PW, R_PW[2,])
    RF_PW = rbind(RF_PW, R_PW[1,])
  }
}

Ridge_PW_new = c()
MLP_PW_new = c()
Lasso_PW_new = c()
ENet_PW_new = c()
RF_PW_new = c()

C = 1
for(l in 1:N_drugs){
  if(l %in% I_zeros){
    
    Ridge_PW_new = rbind(Ridge_PW_new, rep(0,5))
    MLP_PW_new = rbind(MLP_PW_new, rep(0,5))
    Lasso_PW_new = rbind(Lasso_PW_new, rep(0,5))
    ENet_PW_new = rbind(ENet_PW_new, rep(0,5))
    RF_PW_new = rbind(RF_PW_new, rep(0,5))
    
  }else{
    
    Ridge_PW_new = rbind(Ridge_PW_new, Ridge_PW[C,])
    MLP_PW_new = rbind(MLP_PW_new, MLP_PW[C,])
    Lasso_PW_new = rbind(Lasso_PW_new, Lasso_PW[C,])
    ENet_PW_new = rbind(ENet_PW_new, ENet_PW[C,])
    RF_PW_new = rbind(RF_PW_new, RF_PW[C,])
    
    C = C+1
  }
}

rownames(Ridge_PW_new) = colnames(sen)
rownames(MLP_PW_new) = colnames(sen)
rownames(Lasso_PW_new) = colnames(sen)
rownames(ENet_PW_new) = colnames(sen)
rownames(RF_PW_new) = colnames(sen)



# Finding drugs with the same mechanism of actions
Moa = unique(response$moa)
Moa = strsplit(Moa,", ")
Moa = unique(unlist(Moa))

w = strsplit(response$moa, ",")
I_NA = which(is.na(w))
for(r in I_NA){
  w[[r]] = "Unknown"
}

d = list()
l = c()
D = c()
i=0
for(m in Moa){
  print(i)
  i = i+1
  d[m] = list(unique(response[which(sapply(w, FUN=function(x) m %in% x)),12]))
  l = c(l, length(d[[m]]))
  if(length(d[[m]])>5){
    D = c(D, m)
  }
}
#saveRDS(d, "Final_Result/List_of_Drug_MOA.rds")
d = readRDS("Final_Result/List_of_Drug_MOA.rds")
D = c()
for(m in Moa){
  if(length(d[[m]])>5){
  D = c(D, m)
  }
}

# SD = c()
# L = c()
# for(i in 1:length(D)){
#   SD = c(SD, sd(Ridge_PW[d[[D[i]]],1]))
#   L = c(L, length(Ridge_PW[d[[D[i]]],1]))
# }
# hist(SD,100)
# plot(SD,L)
# boxplot(Ridge_Landmark[d[[D[2]]],1])
# 
# max(SD)
#S = which(SD<0.24)
#D_new = D[S]

## obtain data for Landmark
data = list()
M = c()
for(b in 1:length(D)){
  val = Ridge_Landmark[d[[D[b]]],1]

  data[b] = list(val)
  M = c(M, median(val))
}

## obtain data for PW
# data = list()
# M = c()
# for(b in 1:length(D)){
#   val = Ridge_PW_new[d[[D[b]]],1]
#   val_Non_Zeros = val[val!=0]
#   
#   data[b] = list(val_Non_Zeros)
#   M = c(M, median(val_Non_Zeros))
# }

Order = order(M, decreasing = FALSE)
D[Order]

pdf("Figures/FS/Results_whithin/Result5/Boxplot_Moa_Landmark_All_MLs.pdf", height = 10, width = 4)
boxplot(data[Order], horizontal = TRUE)

dev.off()
