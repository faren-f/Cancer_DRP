rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

library(caTools)
source("F10-Ridge.R")
source("F6-ENet.R")
source("F8-MLP.R")
source("F13-Lasso.R")
source("F7-RandomForest.R")

sen = readRDS("Processed_data/S1/sensitivity_matrix_AUC.rds")
TF = read.table("Result_from_Python/TF(gsea2)_PRISM/TF(gsea2)_PRISM_old.csv",
                sep = ",",header = TRUE, row.names = 1)

TF = scale(TF)

#Models = c("RandomForest","ElasticNet", "Lasso","Ridge","MLP")
Models = c("RandomForest", "Ridge")

Mean_Corr_models = c()
STD_Corr_models = c()

Mean_MSE_models = c()
STD_MSE_models = c()

for (M in Models){            # model loop
  model = get(M)
  Mean_Corr = c()
  STD_Corr = c()
  
  Mean_MSE = c()
  STD_MSE = c()
  
  for (i in 86:86){             # drug loop
    print(i)
    Corr = c()
    MSE = c()
    
    for (j in 1:2){           # repeat loop
      print(j)
      X = TF[!is.na(sen[,i]),]
      y = sen[!is.na(sen[,i]),i]
      
      y = scale(y)
      y = y[,1]
      
      sample = sample.split(y, SplitRatio = .8)
      
      Xtrain = subset(X, sample == TRUE)
      Xtest  = subset(X, sample == FALSE)
      ytrain = subset(y, sample == TRUE)
      ytest  = subset(y, sample == FALSE)
      
      
      # Models
      y_pred = model(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
      
      # Evaluation
      corr = cor(ytest,y_pred)
      Corr = c(Corr, corr)
      
      mse =  mean((ytest - y_pred)^2)
      MSE = c(MSE, mse)
    }
    
    Mean_Corr = rbind(Mean_Corr, mean(Corr))
    STD_Corr = rbind(STD_Corr, sd(Corr))
    
    Mean_MSE = rbind(Mean_MSE, mean(MSE))
    STD_MSE = rbind(STD_MSE, sd(MSE))
    
  }
  Mean_Corr_models = cbind(Mean_Corr_models, Mean_Corr)
  STD_Corr_models = cbind(STD_Corr_models, STD_Corr)
  
  Mean_MSE_models = cbind(Mean_MSE_models, Mean_MSE)
  STD_MSE_models = cbind(STD_MSE_models, STD_MSE)
}

# rownames(Mean_Corr_models) = colnames(sen)[1:2]
# colnames(Mean_Corr_models) = Models
# 
# rownames(STD_Corr_models) = colnames(sen)[1:2]
# colnames(STD_Corr_models) = Models
# 
# rownames(Mean_MSE_models) = colnames(sen)[1:2]
# colnames(Mean_MSE_models) = Models
# 
# rownames(STD_MSE_models) = colnames(sen)[1:2]
# colnames(STD_MSE_models) = Models

print(Mean_Corr_models)
system("say Faren Just finished the runing!")

