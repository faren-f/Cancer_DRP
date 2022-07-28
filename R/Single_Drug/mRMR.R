rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Single_Drug/")

## Read data
GE = readRDS("Raw_data/expresion_matrix.rds")
sen = readRDS("Raw_data/sensitivity_matrix.rds")

i=325
Y = sen[1:6,i]
Y = Y[!is.na(sen[,i])]
X = GE[!is.na(sen[,i]),] # remove cell lines that are "NA" For each drug   

#Corr = cor(X,Y)
#high_corr = order(Corr,decreasing = TRUE)
#X = X[,high_corr[1:5]]
X = X[1:6,1:8]

cor_xy = abs(cor(X,Y))
Cor_xx = abs(cor(X,X))

cor_xy_orders = order(cor_xy, decreasing = TRUE)
selected_feat = cor_xy_orders[1]
residual_feat = cor_xy_orders[-1]
scores = max(cor_xy)
alpha = 1

while(length(residual_feat)!=0){
  S = c()
  for(j in residual_feat){
    score_j = cor_xy[j] - alpha*(sum(Cor_xx[j, selected_feat])/length(selected_feat))
    S = c(S, score_j)
  }
  selected_feat = c(selected_feat, residual_feat[which.max(S)])
  residual_feat = residual_feat[-which.max(S)]
  scores = c(scores, max(S))
}
plot(scores)



