rm(list= ls())

library(igraph)
require(caTools)
library(randomForest)


setwd("~/Desktop/Cancer_DRP/R/Single_Drug/Infogenes/")

sen = readRDS("Raw_data/PRISM/sensitivity_matrix.rds")

# STRING
GE = readRDS("Raw_data/PRISM/expresion_matrix_PRISM_STRING.rds")
STRING_edgelist = readRDS("Raw_data/PRISM/ppi_STRING_PRISM.rds")

# Omnipath
#GE = readRDS("Raw_data/PRISM/expresion_matrix_PRISM_Omnipath.rds")
#Omnipath_edgelist = readRDS("Raw_data/PRISM/ppi_Omnipath_PRISM.rds")

# Build graph -------------------------------------------------------------
Interaction_Network = "STRING"

if (Interaction_Network == "STRING"){
  ppi = rbind(STRING_edgelist$gene_symbol1, STRING_edgelist$gene_symbol2)
  
} else if(Interaction_Network == "Ompnipath"){
  ppi = rbind(Omnipath_edgelist$gene_symbol1, Omnipath_edgelist$gene_symbol2)
}

MyGraph = simplify(graph(ppi, directed = FALSE))
my_genes = V(MyGraph)$name


# Remove cell lines that do not have drug response from GE and sen
i = 325                            # drug number
X = GE[!is.na(sen[,i]),]
y = sen[!is.na(sen[,i]),i]

# Normalization
X = scale(X)
y = scale(y)
Rep = 100
mse = rep(0,Rep)
corr = rep(0,Rep)
for (j in 1:Rep){
  print(paste0("Rep is: ",j))
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
  N1 = 100
  N2 = N1+500
  
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
  
  mse[j] = mean((ytest - y_pred)^2)
  corr[j] = cor(ytest, y_pred, method = "pearson")
  
  #print(mse[j])
  print(corr[j])
  
  }

Test_Result = data.frame(mse = mse, corr = corr)
print(apply(Test_Result,2,mean))
#print(apply(Test_Result,2,sd))


