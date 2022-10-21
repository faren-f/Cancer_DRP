rm(list=ls())

source("F10-Ridge.R")
library(parallel)
no_cores = detectCores()
cl = makeCluster(no_cores-2)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
sen = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")
GE = readRDS("Processed_Data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")

source("F14-Feature_Selection.R")
selected_features = c("Landmark_genes")
Omics_List = Feature_Selection(selected_features,GE)
omics = Omics_List[[1]]
index = Omics_List[[2]]

i=24

X = omics[!is.na(sen[,i]),]
y = sen[!is.na(sen[,i]),i]
#hist(y)
#Normalization
X_T = t(X)
X_T = scale(X_T)
X = t(X_T)

clusterExport(cl, c("X","y"))
clusterEvalQ(cl, c(source("F10-Ridge.R")))

RepeatLoop = function(l){
  y_pred_Ridge_all = c()
#for(l in 1:nrow(X)){
  
  Xtrain = X[-l,]
  Xtest = rbind(X[l,],X[l,])
  ytrain = y[-l]

  y_pred_Ridge = Ridge(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
  y_pred_Ridge_all = c(y_pred_Ridge_all, y_pred_Ridge[1])
  return(y_pred_Ridge_all)
}

y_pred_Ridge_all = parLapply(cl, sapply(1:nrow(X), list), RepeatLoop) 

y_pred_Ridge_all = data.frame(y_pred_Ridge_all)

cor(y, y_pred_Ridge_all, method = "pearson")

stopCluster(cl)






badSamples = readRDS("All_Results/XI_Normalized_20%_180&139.rds")
#badSamples1 = readRDS("All_Results/XI_Normalized_20%_D18_90@290.rds")
#badSamples = intersect(badSamples,badSamples1)


color = rep("black", nrow(X))
names(color) = rownames(X)
color[badSamples] = "red"

plot(y, y_pred_Ridge_all, pch = 20, col = color)
# 
# y_new = y[!(rownames(X) %in% intersect(badSamples,rownames(X)))]
# y_pred_Ridge_all_new = y_pred_Ridge_all[!(rownames(X) %in% intersect(badSamples,rownames(X)))]
# 
# cor(y_new, y_pred_Ridge_all_new)
# 
# 


