rm(list = ls())
setwd("~/Desktop/Codes/R/PhD_Project/DRP_PRISM/")

# Library -----------------------------------------------------------------
library(randomForest)
library(ROCR)
library(e1071)

# Read_Data ---------------------------------------------------------------
expr = readRDS("Data/Processed_Data/expresion_matrix.rds")
expr_norm = readRDS("Data/Processed_Data/expresion_normalized_matrix.rds")

sen = readRDS("Data/Processed_Data/sensitivity_matrix.rds")

## Find drug responses for drug i and remove NA samples 
##from sen and expr matrix 

sen_drug_i_with_NA = sen[,20]
indeces_without_NA = which(!is.na(sen_drug_i_with_NA))
sen_drug_i = sen_drug_i_with_NA[indeces_without_NA]
target = sen_drug_i
thr_target = 1.15
hist(target,30)
abline(v =thr_target,col = "red")

target = ifelse(target>thr_target,1,0)
y = as.factor(target+1)

expr_i = expr[indeces_without_NA,]

## To check if there is any cor between response and genes
thr_cor = 0
cor_Xy = apply(expr_i,2,function(x){abs(cor(x,as.numeric(y)))})
hist(cor_Xy,20)

### Feature selection
sd_expr = apply(expr_i, 2, sd)
thr_feat = 1.2
hist(sd_expr,30)
abline(v = thr_feat,col = "red")
ind_high_sd = which(sd_expr>thr_feat)
X = data.frame(expr_i[,ind_high_sd])
#X = expr_i


N_sample = nrow(X)

## Hyperparameters

#kernelParam = 0.08      ## kernelParam gooooood
kernelParam = 0.08       ## 

Rep = 5
Result_AUC_Final = c()
From = 200
To = 200
Step = 1
for (nFeat in seq(From,To,Step)) {
  Result_AUC = c()
  for (rep in 1:Rep){
    y_hat = c()
    
    for (i in 1:N_sample){
      X_test = X[i,]
      y_test = y[i]
      X_train = X[-i,]
      y_train = y[-i]
      
      
      ## Internal Feature selection
      #1)corr with sen
      # nFeat = 150
      FeatureCorrelation = apply(X_train, 2, function(x){abs(cor(x,as.numeric(y_train)))})
      featureSet = order(FeatureCorrelation, decreasing = TRUE)[1:nFeat]
      #hist(featureSet,30)
      
      ##2)p_val
      ## P_val distribution
      #  FeaturePval = apply(X_train, 2, function(x){
      #   a = x[y_train==levels(y_train)[1]]
      #   b = x[y_train==levels(y_train)[2]]
      #   p_val = t.test(a,b)$p.value
      #   return(p_val)
      # })
      #hist(FeaturePval,30)
      #featureSet = order(FeaturePval, decreasing = FALSE)[1:nFeat]
      #Undersample the majarity class
      n1 = table(y_train)[1]
      n2 = table(y_train)[2]
      majorityClass = ifelse(n1>n2,1,2)
      n_resample = abs(n1-n2)
      
      ind_mC = which(y_train == levels(y_train)[majorityClass])
      ind_mC_extra = sample(ind_mC, n_resample, replace = FALSE)
      
      X_train_resampled = X_train[-ind_mC_extra,]
      y_train_resampled = y_train[-ind_mC_extra]
      
      #X_train_resampled = X_train
      #y_train_resampled = y_train
      
      ## Random Forest with internal Feature selection
      model = randomForest(y = as.factor(y_train_resampled),x = X_train_resampled[,featureSet], ntree = 100 ,mtry = 10)
      y_hat_i = predict(model, newdata = X_test[,featureSet])
      
      ## SVM with internal Feature selection
      #model = svm(y = as.factor(y_train_resampled), x = X_train_resampled[,featureSet],kernel = "radial", gamma = kernelParam, scale = T)
      #y_hat_i = predict(model, newdata = rbind(X_test[,featureSet]))
      
      ## SVM without Feature selection
      #model = svm(y = as.factor(y_train_resampled), x = X_train_resampled, kernel = "radial", gamma = kernelParam, scale = F)
      #y_hat_i = predict(model, newdata = rbind(X_test))
      
      y_hat = c(y_hat, as.numeric(y_hat_i))
    }
    
    pred = prediction(y_hat, y)
    auc = performance(pred, measure = "auc")
    
    #auc = performance(pred, "tpr", "fpr")
    
    Result_AUC = c(Result_AUC, as.numeric(auc@y.values))
    #plot(Result_AUC,
    #avg= "threshold", colorize=TRUE, lwd= 3,
    #main= "With ROCR you can produce standard plots\nlike ROC curves ...")
  }
  
  print(Result_AUC)
  print(mean(Result_AUC))
  print(sd(Result_AUC))
  
  Result_AUC_Final = c(Result_AUC_Final, mean(Result_AUC))
  
}

plot(Result_AUC_Final)


