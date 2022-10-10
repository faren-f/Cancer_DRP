rm(list=ls())

library(parallel)
no_cores = detectCores()
cl = makeCluster(no_cores-1)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
sen = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")
GE = readRDS("Processed_Data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")


source("F14-Feature_Selection.R")
selected_features = c("Landmark_genes")
Omics_List = Feature_Selection(selected_features,GE)
omics = Omics_List[[1]]
index = Omics_List[[2]]

N_drug = ncol(sen)
Results = c()
for (i in 5){
  print(paste0("The drug number is: ", as.character(i)))
  
  X = omics[!is.na(sen[,i]),]
  y = sen[!is.na(sen[,i]),i]
  #boxplot(y)
  
  U = 1.5*IQR(y)+median(y)
  L = median(y)-1.5*IQR(y)
  X = X[y<U & y>L,]
  y = y[y<U & y>L]

  
  X = t(scale(t(X)))
  
  
  clusterExport(cl, c("X","y","i","index"))
  clusterEvalQ(cl, c(library(caTools),source("F10-Ridge.R")))
  
  RepLoop = function(j){
    
    sample = sample.split(y, SplitRatio = .9)
    
    Xtrain = subset(X, sample == TRUE)
    Xtest  = subset(X, sample == FALSE)
    ytrain = subset(y, sample == TRUE)
    ytest  = subset(y, sample == FALSE)
    
    
    y_pred_Ridge = Ridge(ytrain = ytrain, Xtrain = Xtrain, Xtest = Xtrain)
    error = (ytrain-y_pred_Ridge)^2
    # mean(error)
    a = hist(error, 10)
    InvCount = 1/a[["counts"]]
    InvCount[is.infinite(InvCount)] = 1
    
    W = rep(0,length(error))
    b = a$breaks
    for(l in 1:length(error)){
      for(k in 2:length(b)){
        if((b[k-1] < error[l]) & (error[l] < b[k])){
          W[l] = InvCount[k-1]
        }
      }
    }
    y_pred_Ridge_after = Ridge(ytrain = ytrain, Xtrain = Xtrain, Xtest = Xtest, weight = W)
    error = (ytest-y_pred_Ridge_after)^2
    mean(error)
    corr_Ridge = cor(ytest,y_pred_Ridge_after,method = "pearson")
    print(corr_Ridge)

    result = data.frame(corr_Ridge = corr_Ridge)
    
    return(result)
  }
  
  N_itration = 100
  result = parLapply(cl, sapply(1:N_itration, list), RepLoop) 
  
  Result = data.frame()
  for (k in 1:N_itration){
    Result = rbind(Result, result[[k]])
  }
  
  Result_mean = apply(Result, 2, mean)
  #Result_sd = apply(Result, 2, sd)
  print(Result_mean)
  #print(Result_sd)
  
  
  Results = rbind(Results, Result_mean)
  
}
stopCluster(cl)



