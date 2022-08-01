rm(list=ls())

require(caTools)
library(randomForest)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

GE = readRDS("Processed_Data/expresion_matrix.rds")
sen = readRDS("Processed_Data/sensitivity_matrix.rds")

# Finding transcription activities using dorothea
data(dorothea_hs, package = "dorothea")

GE = t(GE)              # input: rows are genes and columns are cell lines 
tf_activities <- run_viper(GE, dorothea_hs,
                           options =  list(method = "scale", minsize = 4,
                                           eset.filter = FALSE, cores = 1,
                                           verbose = FALSE))

TF = t(tf_activities)


i = 9                            # drug number
X = TF[!is.na(sen[,i]),]
y = sen[!is.na(sen[,i]),i]

# Normalization
X = scale(X)
y = scale(y)
Rep = 10
MSE = rep(0,Rep)
Corr = rep(0,Rep)
for (j in 1:Rep){
  
  ## Split data into train & test
  sample = sample.split(y, SplitRatio = .8)
  
  Xtrain = subset(X, sample == TRUE)
  Xtest  = subset(X, sample == FALSE)
  ytrain = subset(y, sample == TRUE)
  ytest  = subset(y, sample == FALSE)
  ytrain = as.vector(ytrain)
  ytest = as.vector(ytest)
  ## train model
  RF = randomForest(y = ytrain,x = Xtrain, ntree = 200,mtry = 50)
  y_pred = predict(RF, newdata=Xtest)
  
  MSE[j] = mean((ytest - y_pred)^2)
  Corr[j] = cor(ytest, y_pred, method = "pearson")
  
  #print(MSE[j])
  print(Corr[j])
  
}

Test_Result = data.frame(MSE = MSE, Cor = Corr)
print(apply(Test_Result,2,mean))

