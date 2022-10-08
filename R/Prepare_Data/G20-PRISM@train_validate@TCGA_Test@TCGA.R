rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

source("F7-RandomForest.R")
source("F6-ENet.R")
source("F8-MLP.R")
source("F11-SGL.R")
source("F13-Lasso.R")
source("F18-Combat_Normalization.R")
source("F23-Ridge_TCGA_validation.R")

sen_PRISM = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")
res_TCGA = readRDS("Processed_data/Other/Res_TCGA_24_Drugs.rds")

GE = readRDS("Processed_data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")
GE_TCGA = readRDS("Processed_data/S23/expresion_matrix_TCGA.rds")

# Remove genes whose Q3 is zero
q3_genes = apply(GE_TCGA,2,quantile,prob=0.75)
GE_TCGA = GE_TCGA[,-which(q3_genes==0)]
GE = GE[,-which(q3_genes==0)]

N_drug = ncol(sen_PRISM)
drugs = data.frame(colnames(sen_PRISM))
#Remove:3,10,12,14
#Remove later: 7,11,18,19,24

d = 5
drug = drugs[d,1]

X_PRISM = GE[!is.na(sen_PRISM[,d]),]
y_PRISM = sen_PRISM[!is.na(sen_PRISM[,d]),d]
hist(y_PRISM)
    
X_TCGA = GE_TCGA[!is.na(res_TCGA[,d]),]
y_TCGA = res_TCGA[!is.na(res_TCGA[,d]),d]
    
X_Normalization = Combat_Scale(X_PRISM,X_TCGA)

X_PRISM = X_Normalization[[1]]
X_TCGA = X_Normalization[[2]]
N_genes = ncol(X_PRISM)


source("F15-Feature_Selection_PRISM@TCGA.R")
selected_features = c("Landmark_genes")
Omics_List = Feature_Selection_PRISM_TCGA(selected_features,GE = X_PRISM ,GE_test = X_TCGA)
X_PRISM = Omics_List[[1]]
index = Omics_List[[2]]
X_TCGA = Omics_List[[3]]

length(y_TCGA)
ceiling(0.5*sum(y_TCGA==1)) 
ceiling(0.5*sum(y_TCGA==2)) 

Names_Xval = c(names(sample(y_TCGA[y_TCGA==1],ceiling(0.3*sum(y_TCGA==1)))),
         names(sample(y_TCGA[y_TCGA==2],ceiling(0.3*sum(y_TCGA==2)))))
Xval = X_TCGA[Names_Xval,]
Xtest = X_TCGA[!(rownames(X_TCGA) %in% Names_Xval),]
yval = y_TCGA[Names_Xval]
ytest = y_TCGA[!(rownames(X_TCGA) %in% Names_Xval)]

# Models
#y_pred_Ridge = My_SGL(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest,index = index)
#y_pred_Ridge = RandomForest(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
#y_pred_Ridge = ElasticNet(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
#y_pred_Ridge = Lasso(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
ytrain = y_PRISM
Xtrain = X_PRISM
yval = yval
Xval = Xval
y_pred_Ridge = Ridge(ytrain = y_PRISM ,Xtrain = X_PRISM, yval = yval, Xval = Xval)
#y_pred_Ridge = MLP(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)

# Evaluation
corr_Ridge = cor(y_TCGA,y_pred_Ridge)
#print(corr_Ridge)
#corr_RF = cor(ytest,y_pred_RF)
#corr_ENet = cor(ytest,y_pred_ENet)
#corr_Lasso = cor(ytest,y_pred_Lasso)
#corr_Ridge = cor(ytest , y_pred_Ridge)
#corr_Ridge = cor(ytest , y_pred_Ridge)

#corr_Ridge_PRISM = cor(ytrain_test,y_pred_Ridge_PRISM)
#plot(ytrain_test,y_pred_Ridge_PRISM,xlim = c(0,1.7),ylim = c(0,1.7))


ttest = t.test(y_pred_Ridge[ytest==1], y_pred_Ridge[ytest==2], alternative="greater")$p.value
Ranksum = wilcox.test(y_pred_Ridge[ytest==1], y_pred_Ridge[ytest==2], alternative ="greater")$p.value
#print(Ranksum)
    
plot(ytest,y_pred_Ridge)
boxplot(y_pred_Ridge[ytest==1], y_pred_Ridge[ytest==2])

result = data.frame(corr_Ridge = corr_Ridge, ttest=ttest, 
                    Ranksum = Ranksum)

