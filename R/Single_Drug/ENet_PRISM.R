rm(list = ls())

require(caTools)
library(corrplot)
library(gelnet)
library(RColorBrewer)

setwd("~/Desktop/Cancer_DRP/R/Single_Drug/")

## Read data
GE = readRDS("Raw_data/expresion_matrix.rds")
sen = readRDS("Raw_data/sensitivity_matrix.rds")

## Classifier
Final_Cor = rep(0,20)
mse_test = rep(0,20)

for (j in 1:20){
  i=5
  print(j)
  #Corr = matrix(0,ncol(GE),ncol(sen))
  #for(i in 1:ncol(sen)){
  sen_i = sen[,i]
  sen_i = sen_i[!is.na(sen[,i])]
  
  GE_i = GE[!is.na(sen[,i]),] # remove cell lines that are "NA" For each drug   
  
  ## Corr input & output
  #Corr[,i] = c(cor(GE_i,sen_i))
  #max_cor[i] = max(Corr)  
  
  Corr = cor(GE_i,sen_i)  
  #max(Corr)  
  #hist(abs(Corr)) 
  #}
  ind_Corr = which(abs(Corr)> 0.1)
  sel_GE_i = GE_i[,ind_Corr]
  #expr_each_norm = scale(sel_GE_i)
  #GE_i_norm = (sel_GE_i-min(sel_GE_i))/(max(sel_GE_i)-min(sel_GE_i))
  GE_i_norm = GE_i
  
  #sen/res
  #median_sen_mat = median(sen_i)
  #sen_i_binarized = as.factor(ifelse (sen_i>median_sen_mat,0,1))
  
  # sensitivity normalization
  sen_i_norm = (sen_i-min(sen_i))/(max(sen_i)-min(sen_i))
  
  ## Split data into train & test
  sample = sample.split(sen_i_norm, SplitRatio = .8)
  
  Xtrain = subset(GE_i_norm, sample == TRUE)
  Xtest  = subset(GE_i_norm, sample == FALSE)
  Ytrain = subset(sen_i_norm, sample == TRUE)
  Ytest  = subset(sen_i_norm, sample == FALSE)
  
  
  ## Model training
  n_feat = ncol(GE_i_norm)
  lambda1 = 0.01
  lambda2 = 20
  d = rep(1, n_feat)
  
  GelNet = gelnet(Xtrain, Ytrain, l1 = lambda1, l2 = lambda2, d = d,
                  P = diag(d), m = rep(0,n_feat), max.iter = 50, eps = 1e-05)
  
  ## Test the model
  Beta = GelNet[["w"]]
  Beta0 = GelNet[["b"]]
  y_hat = (Xtest %*%  Beta) + Beta0
  y = Ytest
  
  mse_test[j] = mean((y - y_hat)^2)
  Final_Cor[j] = cor(y,y_hat, method = "pearson")
  
  print(mse_test)
  print(cor(y,y_hat))
}
#}

print(mean(mse_test))
print(mean(Final_Cor))
print(sd(Final_Cor))
