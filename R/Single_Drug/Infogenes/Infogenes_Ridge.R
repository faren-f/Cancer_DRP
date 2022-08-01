rm(list= ls())

library(gelnet)
library(igraph)
require(caTools)

setwd("~/Desktop/Cancer_DRP/R/Single_Drug/Infogenes/")

GE = readRDS("Raw_data/PRISM/expresion_matrix_PRISM_ppi.rds")
sen = readRDS("Raw_data/PRISM/sensitivity_matrix.rds")
ppi_edgelist = readRDS("Raw_data/PRISM/ppi_EdgeList_compelete_PRISM.rds")


# Build graph -------------------------------------------------------------
ppi = rbind(ppi_edgelist$gene_symbol1, ppi_edgelist$gene_symbol2)
MyGraph = simplify(graph(ppi, directed = FALSE))
my_genes = V(MyGraph)$name



# Remove cell lines that do not have drug response from GE and sen
i = 325                            # drug number
X = GE[!is.na(sen[,i]),]
y = sen[!is.na(sen[,i]),i]

# Normalization
X = scale(X)
y = scale(y)
Rep = 10
MSE = rep(0,Rep)
Corr = rep(0,Rep)
ModelParam = data.frame()
for (j in 1:Rep){
  
  ## Split data into train & test
  sample = sample.split(y, SplitRatio = .8)
  
  Xtr_val = subset(X, sample == TRUE)
  Xtest  = subset(X, sample == FALSE)
  ytr_val = subset(y, sample == TRUE)
  ytest  = subset(y, sample == FALSE)
  
  #sample = sample.split(ytr_val, SplitRatio = .8)
  
  #Xtrain = subset(Xtr_val, sample == TRUE)
  #Xval  = subset(Xtr_val, sample == FALSE)
  #ytrain = subset(ytr_val, sample == TRUE)
  #yval  = subset(ytr_val, sample == FALSE)
  
  Xtrain = Xtr_val
  ytrain = ytr_val

  Cor = abs(cor(Xtrain,ytrain))
  Cor_ind = order(Cor,decreasing = TRUE)
  
  # hyperparameters
  N1 = 200
  N2 = N1+100
  
  inf_gene_ind = Cor_ind[1:N1]
  inf_gene = my_genes[inf_gene_ind]
  
  ## Find Diffused genes
  InitialScores = setNames(rep(0,length(my_genes)), my_genes)
  InitialScores[inf_gene]=1
  
  Score_diff = page_rank(MyGraph,directed = FALSE,damping = 0.8,
                personalized = InitialScores)$vector
  Score_diff = Score_diff/max(Score_diff)
  
  S_ind = order(Score_diff, decreasing = TRUE)
  S_ind = S_ind[1:N2]
  good_genes = my_genes[S_ind]
  
  Xtrain = Xtrain[,good_genes]
  Xtest = Xtest[,good_genes]
  
  ## train model
  n_feat = length(good_genes)
  lambda2 = 0.0001
  lambda1 = 0
  d = rep(1, n_feat)
  GelNet = gelnet(Xtrain, ytrain, l1 = lambda1, l2 = lambda2, d = d,
                  P = diag(d), m = rep(0,n_feat), max.iter = 10, eps = 1e-05)
  
  Beta = GelNet[["w"]]
  Beta0 = GelNet[["b"]]
  y_pred = (Xtest %*%  Beta) + Beta0
 
  
  MSE[j] = mean((ytest - y_pred)^2)
  Corr[j] = cor(ytest, y_pred, method = "pearson")
  
  #print(MSE[j])
  print(Corr[j])
  
  
  tempModel = data.frame(N1 = N1, N2 = N2, lambda2 = lambda2, W0 = Beta0)
  Genes = data.frame(rep(0,length(my_genes)))
  row.names(Genes) = my_genes
  
  tempModel = cbind(tempModel,t(Genes))
  tempModel[1,names(Beta)] = Beta
  
  ModelParam = rbind(ModelParam,tempModel)
}

Test_Result = data.frame(MSE = MSE, Cor = Corr)
print(apply(Test_Result,2,mean))
#print(apply(Test_Result,2,sd))

## Save File Path
# if (drugi<10) {
#   DrugNo = paste0("00",as.character(drugi))
# }else if (drugi<100){
#   DrugNo = paste0("0",as.character(drugi))
# }else{
#   DrugNo = as.character(drugi)
# }
# FilePath = paste0("Results_InfoGenes/InfoGenes_",DrugNo,".rds")

#save(Test_Result, ModelParam, file = FilePath)



  