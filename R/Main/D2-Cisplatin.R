rm(list=ls())
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
source("F16-Zscore_Normalization.R")
source("F17-Rank_Normalization.R")
source("F18-Combat_Normalization.R")
source("F7-RandomForest.R")
source("F6-ENet.R")
source("F8-MLP.R")
source("F10-Ridge.R")
source("F11-SGL.R")
source("F13-Lasso.R")
library(caTools)
library(UBL)

cellline_info = read.csv("Raw_data/PRISM/Secondary/secondary-screen-cell-line-info.csv")
Cancer_type = readRDS("Processed_data/S22/Cancer_Types.rds")
sen_PRISM = readRDS("Processed_data/S23/sensitivity_matrix_PRISM_with@TCGA@drugs.rds")
#res_TCGA = readRDS("Processed_data/S24/Drug_response_TCGA_binarized.rds")
res_TCGA = readRDS("Processed_data/S23/Drug_response_matrix_TCGA.rds")
GE = readRDS("Processed_data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")
GE_TCGA = readRDS("Processed_data/S23/expresion_matrix_TCGA.rds")

drugs = data.frame(colnames(sen_PRISM))

# cell_names = rownames(GE)
# cellline_name_tissue = data.frame(cellline_info[cellline_info$row_name %in% cell_names,c(1,4)])
# gastric_celllines = cellline_name_tissue[cellline_name_tissue$primary_tissue=="gastric",1]
# GE = GE[gastric_celllines,]
# sen = sen[gastric_celllines]

source("F15-Feature_Selection_PRISM@TCGA.R")
selected_features = c("Landmark_genes")
Omics_List = Feature_Selection(selected_features,GE = GE ,GE_test = GE_TCGA)
GE = Omics_List[[1]]
index = Omics_List[[2]]
GE_TCGA = Omics_List[[3]]

X = GE[!is.na(sen_PRISM[,53]),]
y = sen_PRISM[!is.na(sen_PRISM[,53]),53]
y= scale(y)
hist(y)
boxplot(y[y<0.5],y[y>0.5])

Results = c()
for(i in 1:5){
  sample = sample.split(y, SplitRatio = .9)
  
  Xtrain = subset(X, sample == TRUE)
  Xtest  = subset(X, sample == FALSE)
  ytrain = subset(y, sample == TRUE)
  ytest  = subset(y, sample == FALSE)
  
  x = data.frame(Xtrain,ytrain)
  balance = SmoteRegress(ytrain~., x, rel = "auto", thr.rel = 0.6, C.perc = "balance",
                          k = 3, repl = FALSE, dist = "Euclidean")
   
  ytrain = balance$ytrain
  #hist(ytrain)
  #hist(ytest)
  Xtrain = balance[,-ncol(balance)]

  
  y_pred_Ridge = Ridge(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
  corr_Ridge = cor(ytest,y_pred_Ridge,method = "pearson")
  corr_Ridge
  #hist(ytrain)
  #hist(ytest)
  #plot(ytest)
  #hist(y_pred_Ridge)
  Results = c(Results,corr_Ridge)
}

print(Results)
mean(Results)
#plot(ytest,y_pred_Ridge)
#error = sqrt((y_pred_Ridge-ytest)^2)
#error = y_pred_Ridge-ytest


# hist(error)
# hist(ytest)
# cor(error,ytest)
# cor(error,y_pred_Ridge)
# hist(y_pred_Ridge)


#########clinical data
GE_TCGA = GE_TCGA[!is.na(res_TCGA[,53]),]
res = res_TCGA[!is.na(res_TCGA[,53]),53]

#GE_TCGA = GE_TCGA[P,]
#res = res[P]



#X_Normalization = Rank(Xtrain,Xtest)
#X_Normalization = Rank(Xtrain,Xtest)
X_Normalization = Combat_Scale(X,GE_TCGA)

X = X_Normalization[[1]]
GE_TCGA = X_Normalization[[2]]

# Models
#y_pred_Ridge = My_SGL(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest,index = index)
#y_pred_RF = RandomForest(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
#y_pred_ENet = ElasticNet(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
#y_pred_Lasso = Lasso(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
y_pred_Ridge = Ridge(ytrain = y ,Xtrain = X, Xtest = GE_TCGA)
#y_pred_Ridge = MLP(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)

# Evaluation
corr_Ridge = cor(res,y_pred_Ridge)
#corr_RF = cor(ytest,y_pred_RF)
#corr_ENet = cor(ytest,y_pred_ENet)
#corr_Lasso = cor(ytest,y_pred_Lasso)
#corr_Ridge = cor(ytest , y_pred_Ridge)
#corr_Ridge = cor(ytest , y_pred_Ridge)

ttest = t.test(y_pred_Ridge[res==1], y_pred_Ridge[res ==4], alternative="greater")$p.value
Ranksum = wilcox.test(y_pred_Ridge[res ==1], y_pred_Ridge[res ==4], alternative ="greater")$p.value


#plot(ytest,y_pred_Ridge)
#boxplot(y_pred_Ridge[ytest==1], y_pred_Ridge[ytest==2])
#plot(ytest,y_pred_SGL)
#corr_MLP = cor(ytest,y_pred_MLP)
result = data.frame(corr_Ridge = corr_Ridge, ttest=ttest, Ranksum = Ranksum)
#corr_SGL = corr_SGL,
#corr_RF = corr_RF,
#corr_ENet = corr_ENet,
#corr_Lasso = corr_Lasso,
#corr_MLP = corr_MLP)
print(result)

e = rownames(res_TCGA)[which(!is.na(res_TCGA[,15]))]
C = Cancer_type[Cancer_type[,2]%in%e,]
C = C[!duplicated(C[,2]),]
P = C[C[,1]=="STAD",2]


boxplot(y_pred_Ridge[res==1], y_pred_Ridge[res==2],
        y_pred_Ridge[res==3], y_pred_Ridge[res==4])







ir <- iris[-c(95:130), ]
hist(x$sen)
mysmote1.iris <- SmoteRegress(Sepal.Width~., ir, dist = "HEOM",
                              C.perc=list(0.5,2.5))

xx <- SmoteRegress(sen~., x, dist = "HEOM",
                              C.perc = "balance", thr.rel = .9)
hist(xx$sen)


smoteBalan.iris <- SmoteRegress(Sepal.Width~., ir, dist = "HEOM",
                                C.perc = "balance")
smoteExtre.iris <- SmoteRegress(Sepal.Width~., ir, dist = "HEOM",
                                C.perc = "extreme")
hist(mysmote2.iris$Sepal.Width)









