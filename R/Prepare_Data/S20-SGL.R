rm(list=ls())
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

source("F4-DoRothEA.R")
source("F9-decoupleR.R")
library(parallel)

no_cores = detectCores()
cl = makeCluster(no_cores-2)
res_drugs = readRDS("Processed_data/S0/Result_All_Drugs.rds")
order_drugs = data.frame(order = order(res_drugs[,3],decreasing = TRUE))

sen = readRDS("Processed_Data/S1/sensitivity_matrix.rds")
GE = readRDS("Processed_Data/S1/expresion_matrix.rds")

# L1000
l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
O1 = GE[,colnames(GE)%in%l1000_genes]

# TF
O2 = DoRothEA(X = GE)

#Tissue types
sample_tissue_types = readRDS("Processed_data/S19/sample_tissue_types.rds")
O3 = sample_tissue_types

# dR
#TF_dR = decoupleR(X = GE, method = "gsva")
O4 = X_TF

# concatenate all omics data
#omics = cbind(O1,O2,O3)
omics = cbind(O1,O2,O4)

#index = c(rep(1,ncol(O1)),rep(2,ncol(O2)),rep(3,ncol(O3)))
index = c(rep(1,ncol(O1)),rep(2,ncol(O2)),rep(2,ncol(O4)))

i = 1399
X = omics[!is.na(sen[,i]),]
y = sen[!is.na(sen[,i]),i]

clusterExport(cl, c("X","y","i","index"))
clusterEvalQ(cl, c(library(caTools),source("F7-RandomForest.R"),
                   source("F11-SGL.R")))

RepLoop = function(j){

  sample = sample.split(y, SplitRatio = .8)
  
  Xtrain = subset(X, sample == TRUE)
  Xtest  = subset(X, sample == FALSE)
  ytrain = subset(y, sample == TRUE)
  ytest  = subset(y, sample == FALSE)
  
  # Normalization-------------------------------------------------------------
  #Xtrain normalization
  Mean_X = apply(Xtrain,2,mean)
  STD_X = apply(Xtrain,2,sd)
  Xtrain = (Xtrain-Mean_X)/STD_X

  # Xtest normalization
  Xtest = (Xtest-Mean_X)/STD_X
  
  #for when we have tissue types
  # Mean_X = apply(Xtrain[,1:2254],2,mean)
  # STD_X = apply(Xtrain[,1:2254],2,sd)
  # Xtrain_1 = (Xtrain[,1:2254]-Mean_X)/STD_X
  # Xtrain = cbind(Xtrain_1,Xtrain[,2255:ncol(Xtrain)])
  # # Xtest normalization
  # Xtest_1 = (Xtest[,1:2254]-Mean_X)/STD_X
  # Xtest = cbind(Xtest_1,Xtest[2255:ncol(Xtest)])
  
  # Ytrain normalization
  # Mean_y = mean(ytrain)
  # STD_y = sd(ytrain)
  # ytrain_norm = (ytrain-Mean_y)/STD_y
  
  # Models
  y_pred_SGL = My_SGL(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest,index = index)
  
  y_pred_RF = RandomForest(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
  
  corr_SGL = cor(ytest,y_pred_SGL)
  corr_RF = cor(ytest,y_pred_RF)
  
  result = data.frame(corr_SGL = corr_SGL, corr_RF = corr_RF)
  
  return(result)
}

N_itration = 6
result = parLapply(cl, sapply(1:N_itration, list), RepLoop) 

Result = data.frame()
for (k in 1:N_itration){
  Result = rbind(Result, result[[k]])
}

r = c(mean(Result$corr_SGL), mean(Result$corr_RF))
print(r)

stopCluster(cl)

