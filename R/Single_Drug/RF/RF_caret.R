rm(list = ls())
library(randomForest)
require(caTools)
library(caret)
setwd("~/Desktop/Cancer_DRP/R/Single_Drug/RF/")

## Read data
GE = readRDS("Raw_data/expresion_matrix.rds")
TF = read.table("Raw_data/TF_gsea_1.csv",header = TRUE,sep = ",",row.names=1)

GE = TF

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
Rep = 1
corr = rep(0,Rep)
mse = rep(0,Rep)

for (j in 1:Rep){
  print(j)
  #i = 540
  i = 1432
  #Corr = matrix(0,ncol(GE),ncol(sen))
  #max_cor = rep(0,ncol(sen))
  #for(i in 1:ncol(sen)){
  y = sen[,i]
  y = y[!is.na(sen[,i])]
  X = GE[!is.na(sen[,i]),] # remove cell lines that are "NA" For each drug   
  Corr = cor(X,y)
  high_corr = order(Corr,decreasing = TRUE)
  #X = X[,high_corr[1:200]]
  
  
  ## Corr input & output
  #Corr[,i] = c(cor(X,y))
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
  #median_sen_mat = median(y)
  #y_binarized = as.factor(ifelse (y>median_sen_mat,1,2))
  
  # sensitivity normalization
  #y_norm = (y-min(y))/(max(y)-min(y))
  y_norm = scale(y)
  ## Split data into train & test
  sample = sample.split(y_norm, SplitRatio = .8)
  
  Xtrain = subset(X_norm, sample == TRUE)
  Xtest  = subset(X_norm, sample == FALSE)
  ytrain = subset(y_norm, sample == TRUE)
  ytest  = subset(y_norm, sample == FALSE)
  
  ytrain = as.vector(ytrain)
  ytest = as.vector(ytest)
  
  
  train_data = cbind(Xtrain,ytrain)
  
  control = trainControl(method = "cv",
                          number = 10,
                          search = "grid",
                         allowParallel = TRUE)
  
  tuneGrid = expand.grid(.mtry = seq(50, 200, by = 5))

  maxtrees = list()
  for (ntree in c(50,70,100,150,200,250)) {
    print(paste0("ntree is: ", ntree))
    model <- caret::train(ytrain~., 
                         data = train_data,
                         method = "rf",
                         tuneGrid = tuneGrid,
                         trControl = control,
                         ntree = ntree)
    key <- toString(ntree)
    maxtrees[[key]] = model
  }
  results_tree <- resamples(maxtrees)
  summary(results_tree)

  # model = caret::train(ytrain ~ ., data = train_data,
  #                method = "rf",
  #                ntree = 100,
  #                trControl = control,
  #                tuneGrid = NULL)
  # 
  y_pred = predict(model, Xtest)
  

  #AC = sum((as.numeric(prediction)==1 & as.numeric(ytest)==1) |
  #(as.numeric(prediction)==2 & as.numeric(ytest)==2))
  #ACC = AC/length(ytest)*100
  
  corr[j] = cor(ytest,y_pred)
  mse[j] = mean((ytest-y_pred)^2)
  
  print(corr)
  #print(mse)
}
#}
#print(corr)
#print(mse)
print(mean(corr))
#print(sd(corr))
#print(mean(mse))
#plot(corr,mse)



