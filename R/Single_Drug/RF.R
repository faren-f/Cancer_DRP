rm(list = ls())

library(randomForest)
require(caTools)

setwd("~/Desktop/Cancer_DRP/R/Single_Drug/")

## Read data
GE = readRDS("Data/Processed_Data/expresion_matrix.rds")
sen = readRDS("Data/Processed_Data/sensitivity_matrix.rds")

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

for (j in 1:20){
  print(j)
  i=1429
  #Corr = matrix(0,ncol(GE),ncol(sen))
  #max_cor = rep(0,ncol(sen))
  #for(i in 1:ncol(sen)){
  sen_i = sen[,i]
  sen_i = sen_i[!is.na(sen[,i])]
  expr_each_drug = GE[!is.na(sen[,i]),] # remove cell lines that are "NA" For each drug   
  Corr = cor(expr_each_drug,sen_i)
  
  ## Corr input & output
    #Corr[,i] = c(cor(expr_each_drug,sen_i))
    #max_cor[i] = max(Corr)  
    #hist(abs(Corr))  
  #}
  
  # e = abs(Corr)
  # good_drug = apply(e,2,sum)
  # d = sort(good_drug, decreasing = TRUE) 
  # d = data.frame(d)  
  # 
    
  #ind_Corr = which(abs(Corr)> 0.2)
  #sel_GE_i = GE_i[,ind_Corr]
  #GE_i_norm = scale(GE_i)
  GE_i_norm = GE_i
  #GE_i_norm = (sel_GE_i-min(sel_GE_i))/(max(sel_GE_i)-min(sel_GE_i))
  
  
  #sen/res
  #median_sen_mat = median(sen_i)
  #sen_i_binarized = as.factor(ifelse (sen_i>median_sen_mat,1,2))
  
  # sensitivity normalization
  sen_i_norm = (sen_i-min(sen_i))/(max(sen_i)-min(sen_i))
  
  ## Split data into train & test
  sample = sample.split(sen_i_norm, SplitRatio = .8)
  
  Xtrain = subset(GE_i_norm, sample == TRUE)
  Xtest  = subset(GE_i_norm, sample == FALSE)
  Ytrain = subset(sen_i_norm, sample == TRUE)
  Ytest  = subset(sen_i_norm, sample == FALSE)
  
  RF = randomForest(y = Ytrain,x = Xtrain, ntree = 200,mtry = 100)
  y_hat = predict(RF, newdata=Xtest)
  #AC = sum((as.numeric(prediction)==1 & as.numeric(Ytest)==1) |
             #(as.numeric(prediction)==2 & as.numeric(Ytest)==2))
  #ACC = AC/length(Ytest)*100
  
  Final_Cor[j] = cor(Ytest,y_hat)
  MSE[j] = mean((Ytest-y_hat)^2)
  
  print(Final_Cor)[j]
  print(MSE)[j]
}
#}
print(Final_Cor)
print(MSE)
print(mean(Final_Cor))
print(sd(Final_Cor))
#print(dim(expr_each_drug_norm))
print(mean(MSE))
#plot(Final_Cor,MSE)



