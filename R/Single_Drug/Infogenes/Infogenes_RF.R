rm(list= ls())

library(igraph)
require(caTools)
library(randomForest)


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
  #yval = as.vector(yval)
  
  Xtrain = Xtr_val
  ytrain = ytr_val
  
  ytrain = as.vector(ytrain)
  ytest = as.vector(ytest)

  Cor = abs(cor(Xtrain,ytrain))
  Cor_ind = order(Cor,decreasing = TRUE)
  
  # hyperparameters
  N1 = 300
  N2 = N1+200
  
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
  
  RF = randomForest(y = ytrain,x = Xtrain, ntree = 200,mtry = 100)
  y_pred = predict(RF, newdata=Xtest)
  
  
  
  
  MSE[j] = mean((ytest - y_pred)^2)
  Corr[j] = cor(ytest, y_pred, method = "pearson")
  
  #print(MSE[j])
  print(Corr[j])
  
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



