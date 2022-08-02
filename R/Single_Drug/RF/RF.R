rm(list = ls())
library(randomForest)
require(caTools)

setwd("~/Desktop/Cancer_DRP/R/Single_Drug/RF/")

## Read data
GE = readRDS("Raw_data/expresion_matrix.rds")
sen = readRDS("Raw_data/sensitivity_matrix.rds")

## Pre Processing
#Remove Genes with Low Median
#median_genes = apply(GE,2,median)
#names(median_genes) = colnames(GE)
#median_genes = sort(median_genes)
#hist(median_genes,100)
#abline(v=3.5,col="red")
#median_genes = median_genes[median_genes>3.5]
#GE = GE[,names(median_genes)]

# Remove Genes with Low STD
# std_genes = apply(GE,2,sd)
# names(std_genes) = colnames(GE)
# std_genes = sort(std_genes)
# hist(std_genes,100)
# abline(v=.4,col="red")
# std_genes = std_genes[std_genes>0.4]
# GE = GE[,names(std_genes)]


## Classifier
Final_Cor = rep(0,20)
MSE = rep(0,20)

for (j in 1:5){
  print(j)
  #i = 541
  i = 265
  #i = 325
  #Corr = matrix(0,ncol(GE),ncol(sen))
  #max_cor = rep(0,ncol(sen))
  #for(i in 1:ncol(sen)){
  Y = sen[,i]
  Y = Y[!is.na(sen[,i])]
  X = GE[!is.na(sen[,i]),] # remove cell lines that are "NA" For each drug   
  Corr = cor(X,Y)
  high_corr = order(Corr,decreasing = TRUE)
  X = X[,high_corr[1:200]]

  
    ## Corr input & output
    #Corr[,i] = c(cor(X,Y))
    #max_cor[i] = max(Corr)  
    #hist(abs(Corr))  
  #}
  
  # e = abs(Corr)
  # good_drug = apply(e,2,sum)
  # d = sort(good_drug, decreasing = TRUE) 
  # d = data.frame(d)  
  # 
    
  #ind_Corr = which(abs(Corr)> 0.2)
  #sel_X = X[,ind_Corr]
  X_norm = scale(X)
  #X_norm = X
  #X_norm = (sel_X-min(sel_X))/(max(sel_X)-min(sel_X))
  
  
  #sen/res
  #median_sen_mat = median(Y)
  #Y_binarized = as.factor(ifelse (Y>median_sen_mat,1,2))
  
  # sensitivity normalization
  #Y_norm = (Y-min(Y))/(max(Y)-min(Y))
  Y_norm = scale(Y)
  ## Split data into train & test
  sample = sample.split(Y_norm, SplitRatio = .8)
  
  Xtrain = subset(X_norm, sample == TRUE)
  Xtest  = subset(X_norm, sample == FALSE)
  Ytrain = subset(Y_norm, sample == TRUE)
  Ytest  = subset(Y_norm, sample == FALSE)
  
  Ytrain = as.vector(Ytrain)
  Ytest = as.vector(Ytest)
  
  RF = randomForest(y = Ytrain,x = Xtrain, ntree = 200,mtry = 100)
  y_pred = predict(RF, newdata=Xtest)
  #AC = sum((as.numeric(prediction)==1 & as.numeric(Ytest)==1) |
             #(as.numeric(prediction)==2 & as.numeric(Ytest)==2))
  #ACC = AC/length(Ytest)*100
  
  Final_Cor[j] = cor(Ytest,y_pred)
  MSE[j] = mean((Ytest-y_pred)^2)
  
  print(Final_Cor)[j]
  print(MSE)[j]
}
#}
print(Final_Cor)
#print(MSE)
print(mean(Final_Cor))
#print(sd(Final_Cor))
#print(dim(X_norm))
#print(mean(MSE))
#plot(Final_Cor,MSE)



