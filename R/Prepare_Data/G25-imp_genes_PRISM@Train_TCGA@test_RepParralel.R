rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
source("F10-Ridge.R")
source("F18-Combat_Normalization.R")

sen_PRISM = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")
res_TCGA = readRDS("Processed_data/Other/Res_TCGA_24_Drugs.rds")

GE_PRISM = readRDS("Processed_data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")
GE_TCGA = readRDS("Processed_data/S23/expresion_matrix_TCGA.rds")

# Remove genes whose Q3 is zero
q3_genes = apply(GE_TCGA,2,quantile,prob=0.75)
GE_TCGA = GE_TCGA[,-which(q3_genes==0)]
GE_PRISM = GE_PRISM[,-which(q3_genes==0)]

i=4
Xtrain = GE_PRISM[!is.na(sen_PRISM[,i]),]
ytrain = sen_PRISM[!is.na(sen_PRISM[,i]),i]

Xtest = GE_TCGA[!is.na(res_TCGA[,i]),]
ytest = res_TCGA[!is.na(res_TCGA[,i]),i]

X_Normalization = Combat_Scale(Xtrain,Xtest)
Xtrain = X_Normalization[[1]]
Xtest = X_Normalization[[2]]

# source("F15-Feature_Selection_PRISM@TCGA.R")
# selected_features = c("Landmark_genes")
# Omics_List = Feature_Selection_PRISM_TCGA(selected_features,Xtrain = Xtrain ,Xtest = Xtest)
# Xtrain = Omics_List[[1]]
# index = Omics_List[[2]]
# Xtest = Omics_List[[3]]


train_data = cbind(Xtrain,ytrain)
control = trainControl(method = "repeatedcv",
                       number = 5,
                       repeats = 5,
                       verboseIter = FALSE)

tune = expand.grid(alpha = 0,lambda = seq(0.01,5,by = 0.1))

model = caret::train(ytrain ~., data = train_data,
                     method = "glmnet",
                     weights = NULL,
                     metric="RMSE",
                     allowParallel = TRUE,
                     tuneGrid = tune,
                     trControl = control)

y_pred = predict(model,Xtest)
Beta = as.matrix(coef(model$finalModel, model$bestTune$lambda))
Order = order(abs(Beta), decreasing = TRUE)

Order_Beta = rep(0,length(Order))
for(j in 1:length(Order)){
  Order_Beta[Order[j]]=j
}

corr = cor(ytest,y_pred)
ttest = t.test(y_pred[ytest==1], y_pred[ytest==2], alternative="greater")$p.value
Ranksum = wilcox.test(y_pred[ytest==1], y_pred[ytest==2], alternative ="greater")$p.value

Result = data.frame(corr = corr, ttest=ttest, Ranksum = Ranksum)
print(Result)

saveRDS(Result,"Final_Result/imp_genes_PRISM&TCGA/Docetaxel/Result_PRISM.rds")
saveRDS(Order_Beta,"Final_Result/imp_genes_PRISM&TCGA/Docetaxel/Order_Beta_PRISM.rds")

