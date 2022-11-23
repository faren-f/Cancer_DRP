rm(list=ls())

library(mnda)
library(keras)
library(aggregation)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
data = readRDS("MNDA/data_lung_tamoxifen.rds")
X = data[[1]]
y = data[[2]]

l_res = sum(y=="res")
l_nonres = sum(y=="non_res")
Cor_res = abs(cor(X[1:l_res,]))
Cor_nonres = abs(cor(X[(l_res+1):(l_res+l_nonres),]))

adj_list = list(Cor_res, Cor_nonres)

graph_data = as.mnda.graph(adj_list, outcome = c("res","non_res"))
emb_list = mnda_embedding_2layer(graph_data, train.rep = 10, walk.rep = 100, n.steps = 2)
Results = mnda_node_detection_2layer(emb_list, p.adjust.method = "none")

