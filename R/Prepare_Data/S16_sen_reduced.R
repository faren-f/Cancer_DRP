rm(list=ls())
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

sen = readRDS("Processed_Data/S1/sensitivity_matrix_AUC.rds")
GE = readRDS("Processed_Data/S1/expresion_matrix.rds")
gene_name = colnames(GE)
corr = matrix(0,ncol(sen),ncol(GE))
for (i in 1:ncol(sen)){
  X = GE[!is.na(sen[,i]),]
  y = sen[!is.na(sen[,i]),i]
  if(nrow(X)>200)
    corr[i,] = abs(cor(X,y))
  else
    corr[i,] = 0
}

corr[is.na(corr)] = 0

saveRDS(corr,"Processed_data/Other/corr_all_drugs_with_sen.rds")


Mean = apply(corr,1,mean)
hist(Mean)
Median = apply(corr,1,median)
hist(Mean)
Q3 = apply(corr,1, function(x){return(quantile(x,0.75))}) 
hist(Q3)
abline(v=.1,col="red")
length(which(Q3>0.1))
Max = apply(corr,1,max)
hist(Max)
abline(v=.2,col="red")
length(which(Max>0.25))
sen_reduced = sen[,which(Max>0.25)]
saveRDS(sen_reduced,"Processed_data/Other/sen_reduced.rds")


res_drugs = readRDS("Processed_data/Other/Result_All_Drugs.rds")
order_drugs = data.frame(order = order(res_drugs[,3],decreasing = TRUE))
#drugs = data.frame(colnames(sen))

which(Max>0.8) 
corr_398 = corr[398,]
which(corr_398>0.8)  
hist(corr_398)  
  
#----------
X = GE[!is.na(sen[,i]),]
y = sen[!is.na(sen[,i]),i]
c = abs(cor(X,y))
hist(c)

X_r = X[1:200,]
y_r = y[1:200]
cc = abs(cor(X_r,y_r))
hist(cc)




  