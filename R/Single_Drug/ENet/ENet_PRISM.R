rm(list = ls())

require(caTools)
library(corrplot)
library(gelnet)
library(RColorBrewer)

setwd("~/Desktop/Cancer_DRP/R/Single_Drug/ENet/")

## Read data
GE = readRDS("Raw_data/expresion_matrix.rds")
sen = readRDS("Raw_data/sensitivity_matrix.rds")

i=325
X = GE[!is.na(sen[,i]),]           # remove cell lines that are "NA" For each drug   
y = sen[!is.na(sen[,i]),i]

# Corr = cor(X,y)
# order_Corr = order(Corr, decreasing = TRUE)
# X = X[,order_Corr[1:200]]

Rep = 2
corr = rep(0,Rep)
mse = rep(0,Rep)
for (j in 1:Rep){
  print(j)
  
  ## Split data into train & test
  sample = sample.split(y, SplitRatio = .8)
  
  Xtrain = subset(X, sample == TRUE)
  Xtest  = subset(X, sample == FALSE)
  ytrain = subset(y, sample == TRUE)
  ytest  = subset(y, sample == FALSE)
    
  # Normalization
  # Xtrain normalization
  Mean_X = apply(Xtrain,2,mean)
  STD_X = apply(Xtrain,2,sd)
  Xtrain = (Xtrain-Mean_X)/STD_X
  
  # Xtest normalization
  Xtest = (Xtest-Mean_X)/STD_X 
  
  #Ytrain normalization
  #Ytrain = (Ytrain-min(Ytrain))/(max(Ytrain)-min(Ytrain))
  Mean_y = mean(ytrain)
  STD_y = sd(ytrain)
  ytrain_norm = (ytrain-Mean_y)/STD_y
  
  ## Model training
  n_feat = ncol(X)
  lambda1 = 0.01
  lambda2 = 20
  d = rep(1, n_feat)
  
  GelNet = gelnet(Xtrain, ytrain, l1 = lambda1, l2 = lambda2, d = d,
                  P = diag(d), m = rep(0,n_feat), max.iter = 50, eps = 1e-05)
  
  ## Test the model
  Beta = GelNet[["w"]]
  Beta0 = GelNet[["b"]]
  y_pred = (Xtest %*%  Beta) + Beta0
  
  # y_pred re-normalization
  y_pred = (y_pred*STD_y)+Mean_y
  
  mse[j] = mean((ytest - y_pred)^2)
  corr[j] = cor(ytest,y_pred, method = "pearson")
  
  print(mse)
  print(corr)
}

#print(mean(mse))
print(mean(corr))
#print(sd(corr))
