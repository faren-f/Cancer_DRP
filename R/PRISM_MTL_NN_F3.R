rm(list = ls())

library(GEOquery)
require(caTools)
library(keras)
library(tensorflow)
library(corrplot)
library(ggvis)

setwd("~/Desktop/Codes/R/PhD_Project/DRP_PRISM/")

#gseID = "GSE36139"
#data = getGEO(gseID)[[1]]
#saveRDS(data, paste0("data/",gseID,".rds"))

## Read data

# Read_Data ---------------------------------------------------------------
GE = readRDS("Data/Processed_Data/expresion_matrix.rds")
drug_sensitivity = readRDS("Data/Processed_Data/sensitivity_matrix.rds")

# Drug features
Fingerprints = readRDS("Data/Processed_Data/Fingerprints.rds")  

FP = matrix(0, length(Fingerprints), 1024)
for (i in 1:length(Fingerprints)){
  FP_bits_on = Fingerprints[[i]]
  FP[i,FP_bits_on@bits] = 1
}
## Pre Processing
# Remove Genes with Low Median
#median_genes = apply(GE,2,median)
#names(median_genes) = colnames(GE)
#median_genes = sort(median_genes)
#hist(median_genes,100)
#abline(v=3.5,col="red")
#median_genes = median_genes[median_genes>0]
#GE = GE[,names(median_genes)]

# Remove Genes with Low STD
#std_genes = apply(GE,2,sd)
#names(std_genes) = colnames(GE)
#std_genes = sort(std_genes)
#hist(std_genes,100)
#abline(v=.4,col="red")
#std_genes = std_genes[std_genes>.1]
#exprData = exprData[,names(std_genes)]

# Clustering the drugs
# Dis_Cor_sen = as.dist(1-abs(cor(sen_mat)))
# Tree = hclust(Dis_Cor_sen,method = "ward.D2")
# plot(Tree)

## Selection of Some Expressions based on (Corr input & output)
#Corr = cor(exprData,sen_mat) 
#Max_Corr = apply(abs(Corr),1,max)
#ind_Corr = which(Max_Corr> 0.1)
#selected_expr = exprData[,ind_Corr]
#selected_expr = exprData

## Normalization of Expressions
#GE_norm = apply(GE,2, function(x){return((x-min(x))/(max(x)-min(x)))})
GE_norm = scale(GE)

# Normalization of Sensitivities
#sen_norm = apply(drug_sensitivity,2, function(x){return((x-min(x))/(max(x)-min(x)))})
sen_norm = scale(drug_sensitivity)

## Classifier
drug_No = ncol(sen_norm)
N = 1
Cor_drugs = matrix(0,drug_No,N)
MSE_drugs = matrix(0,drug_No,N)
for (j in 1:N){
  print(j)
  
  sample = sample.split(sen_norm[,1], SplitRatio = .8)
  
  Xtrain = subset(GE_norm, sample == TRUE)
  Xtest  = subset(GE_norm, sample == FALSE)
  Ytrain = subset(sen_norm, sample == TRUE)
  Ytest  = subset(sen_norm, sample == FALSE)
  
  Ytrain = as.matrix(Ytrain)
  
  ## NA removal
  mean_sen = rep(0,ncol(Ytrain))
  for (k in 1:ncol(Ytrain)){
    
    col = Ytrain[!is.na(Ytrain[,k]),k]
    mean_sen[k] = mean(col)
    Ytrain[is.na(Ytrain[,k]),k] = mean_sen[k]
  }
  
  ###Model
  # Initialize a sequential model
  model <- keras_model_sequential()
  
  # Add layers to the model
  model %>% 
    layer_dense(units = 250, activation = 'relu', input_shape = ncol(Xtrain),
                kernel_regularizer = regularizer_l2(0.001)) %>%  
    #layer_dropout(rate = 0.4) %>% 
    #kernel_initializer = "normal"  
    #bias_initializer = "Zeros") 
    
    layer_dense(units = 100, activation = 'relu') %>% 
    #layer_dropout(rate = 0.3) %>% 
    
    #layer_dense(units = 50, activation = 'relu') %>% 
    #layer_dropout(rate = 0.2) %>% 
    
    
    #layer_dense(units = 10, activation = 'relu') %>%
    layer_dense(units = drug_No, activation = 'linear')
  
  # Compile the model
  model %>% compile(
    loss = 'mse',
    optimizer = "adam")
  
  # Fit the model 
  model %>% fit(
    Xtrain, 
    Ytrain, 
    epochs = 15, 
    batch_size = 5, 
    validation_split = 0.2)
  
  
  # Evaluate on test data and labels
  Y_hat <- model %>% predict(Xtest)
  for (k in 1: ncol(Ytest)){
    
    Ytest2 = Ytest[!is.na(Ytest[,k]),k] 
    Y_hat2 = Y_hat[!is.na(Ytest[,k]),k]
    Cor_drugs[k,j] = cor(Ytest2,Y_hat2)
    MSE_drugs[k,j] = mean((Ytest2-Y_hat2)^2)
  }
}

Final_Cor = apply(Cor_drugs, 1, mean)
Final_MSE = apply(MSE_drugs, 1, mean)

print(Final_MSE)
print(Final_Cor)

