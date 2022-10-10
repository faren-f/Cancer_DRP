rm(list=ls())

library(parallel)
no_cores = detectCores()
cl = makeCluster(no_cores-1)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
source("F18-Combat_Normalization.R")

sen_PRISM = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")
res_TCGA = readRDS("Processed_data/Other/Res_TCGA_24_Drugs.rds")

GE = readRDS("Processed_data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")
GE_TCGA = readRDS("Processed_data/S23/expresion_matrix_TCGA.rds")

# Remove genes whose Q3 is zero
q3_genes = apply(GE_TCGA,2,quantile,prob=0.75)
GE_TCGA = GE_TCGA[,-which(q3_genes==0)]
GE = GE[,-which(q3_genes==0)]


d = 4
X_PRISM = GE[!is.na(sen_PRISM[,d]),]

X_TCGA = GE_TCGA[!is.na(res_TCGA[,d]),]
y_TCGA = res_TCGA[!is.na(res_TCGA[,d]),d]

X_Normalization = Combat_Scale(X_PRISM,X_TCGA)

X_PRISM = X_Normalization[[1]]
X_TCGA = X_Normalization[[2]]


# source("F15-Feature_Selection_PRISM@TCGA.R")
# selected_features = c("Landmark_genes")
# Omics_List = Feature_Selection_PRISM_TCGA(selected_features, Xtrain = X_PRISM, Xtest = X_TCGA)
# X_PRISM = Omics_List[[1]]
# index = Omics_List[[2]]
# X_TCGA = Omics_List[[3]]

length(y_TCGA)
sum(y_TCGA==1)
sum(y_TCGA==2)
X = X_TCGA
y = as.factor(y_TCGA)


clusterExport(cl, c("X","y"))
clusterEvalQ(cl, c(library(ROCR), source("F24-LogisticRegresion.R")))


RepLoop = function(r){
  
  y_pred_logit_all = c()
  Order_Beta_all = c()
  
  for (i in 1:nrow(X)){
    Xtest = X[i,]
    ytest = y[i]
    Xtrain = X[-i,]
    ytrain = y[-i]
    
    #Undersample the majarity class
    n1 = table(ytrain)[1]
    n2 = table(ytrain)[2]
    majorityClass = ifelse(n1>n2,1,2)
    n_resample = abs(n1-n2)
    
    if(n_resample != 0){
      ind_mC = which(ytrain == levels(ytrain)[majorityClass])
      ind_mC_extra = sample(ind_mC, n_resample, replace = FALSE)
      
      Xtrain_resampled = Xtrain[-ind_mC_extra,]
      ytrain_resampled = ytrain[-ind_mC_extra]
      
    }else{
      Xtrain_resampled = Xtrain
      ytrain_resampled = ytrain
    }
    train_data = cbind(Xtrain_resampled,ytrain_resampled)
    
    control = trainControl(method = 'repeatedcv',
                           number = 5,
                           repeats =  5)
    
    model = caret::train(as.factor(ytrain_resampled) ~., data = train_data,
                         method = 'glmnet',
                         trControl = control,
                         allowParallel = TRUE,
                         family = 'binomial')
    
    y_pred = predict(model,rbind(Xtest))
    Beta = as.matrix(coef(model$finalModel, model$bestTune$lambda))
    Order = order(abs(Beta), decreasing = TRUE)
    
    Order_Beta = rep(0,length(Order))
    for(j in 1:length(Order)){
      Order_Beta[Order[j]]=j
    }
    
    Order_Beta_all = cbind(Order_Beta_all,Order_Beta)
    y_pred_logit_all = c(y_pred_logit_all, y_pred)
    
  }
  
  acc = as.numeric(y_pred_logit_all) - as.numeric(y)
  acc = sum(acc == 0)/length(acc)
  
  pred = prediction(y_pred_logit_all, y)
  auc = performance(pred, measure = "auc")
  auc = as.numeric(auc@y.values)
  
  pval = chisq.test(table(data.frame(cbind(y,y_pred_logit_all))))$p.value
  

  result = cbind(acc,auc,pval)
  results = list(result,Order_Beta_all)
  return(results)
}

N_itration = 14
results = parLapply(cl, sapply(1:N_itration, list), RepLoop) 


Result = c()
Order_Beta_all = c()
for (k in 1:N_itration){
  Result = rbind(Result, results[[k]][[1]])
  Order_Beta_all = cbind(Order_Beta_all, results[[k]][[2]])
}

stopCluster(cl)

RankSum_Beta = apply(Order_Beta_all,1,sum)
plot(sort(RankSum_Beta), pch=20)
order(RankSum_Beta)[1:20]

#saveRDS(Result,"Final_Result/imp_genes_PRISM&TCGA/Docetaxel/Result_TCGA.rds")
#saveRDS(Order_Beta_all,"Final_Result/imp_genes_PRISM&TCGA/Docetaxel/Order_Beta_TCGA.rds")

saveRDS(Result,"Final_Result/imp_genes_PRISM&TCGA/Docetaxel/Result_TCGA_WholeGenes.rds")
saveRDS(Order_Beta_all,"Final_Result/imp_genes_PRISM&TCGA/Docetaxel/Order_Beta_TCGA_WholeGenes.rds")

