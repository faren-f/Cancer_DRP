rm(list=ls())

###### Wrapper method
#### using particle swarm optimisation we find the best number of informative genes
####(Correlation or Mutual information) and their neighbours in string with the methods
####of page rank or diffuse and also the best parameters for ENET simultaneously


setwd("~/Desktop/Cancer_DRP/R/Single_Drug/Infogenes/")

library(gelnet)
library(igraph)
library(pso)

##### Load Data
load("Raw_data/GDSC_Expr_IC50_id.rds")
ppi = readRDS("Raw_data/ppi/ppi_gene_900.rds")

Expr = data.frame(Expr)

sen = IC50
n_drug = ncol(IC50)

### Build Graph
ppi = rbind(as.character(ppi[,1]) , as.character(ppi[,2]))
MyGraph = simplify(graph(ppi, directed = FALSE))

### Detect non-intersected genes
Intersect = intersect(colnames(Expr), V(MyGraph)$name)

V_delete = V(MyGraph)$name [! V(MyGraph)$name %in% Intersect] # the nodes that do not exists in intersect
MyGraph = delete_vertices(MyGraph, V_delete)

Expr = as.matrix(Expr[,V(MyGraph)$name])

my_genes = V(MyGraph)$name

## Drug selection loop       
drugi = 1
for (drugi in 1:1) {
  
  ## Find cell lines that have gone under the drug i test
  NA_indx = which (sen[,drugi] %in% NA )
  Expr_i = as.matrix(Expr[-NA_indx,])
  sen_i = as.matrix(sen[-NA_indx,drugi])
  
  ### Normalization
  Expr_i_norm = scale(Expr_i,center = apply(Expr_i,2,mean),
                      scale = apply(Expr_i,2,sd))  
  
  sen_i_norm = scale(sen_i,center = apply(sen_i,2,mean,na.rm = TRUE), 
                     scale = apply(sen_i,2,sd,na.rm = TRUE))
  
  ## Correlation between genes and IC50s
  Corr = abs(apply(Expr_i_norm,2, function(X){return(cor(X,sen_i_norm))}))
  Cor_ind = order(Corr,decreasing = TRUE)
  #Cor_ind = sample(ncol(Expr_i_norm))            # RANDOM TEST
  
  
  ## Change names to X and ic50
  X = data.frame(Expr_i_norm)
  ic50 = sen_i_norm
  
  ## Cross validation loop
  Rep = 1                      # Random sub-sampling repeat
  mse_test = rep(0,Rep)
  corr_test = rep(0,Rep)
  ModelParam = data.frame()
  j=1
  for (j in 1:Rep) {
    
    print(drugi)
    print(j)
    
    ## Train, Validation, Test 
    N = dim(X)[1]                         # Number of cell lines
    a = .6
    b = .8
    Index = sample(1:N)
    Ind_tr = Index[1:floor(a*N)]
    Ind_val = Index[(floor(a*N)+1):floor(b*N)]
    Ind_te = Index[(floor(b*N)+1):N]
    
    X_tr = as.matrix(X[Ind_tr,])          # x train
    X_val = as.matrix(X[Ind_val,])        # x validation
    X_te = as.matrix(X[Ind_te,])          # x test
    
    ic50_tr = ic50[Ind_tr]                # train
    ic50_val = ic50[Ind_val]              # validation
    ic50_te = ic50[Ind_te]                # test
    
    #### Train model using cross validation
    
    ### PSO
    FitnessFunc = function(Param) {
      
      N1 = floor(Param[1])
      N2 = floor(Param[2])+N1
      lambda2 = exp(Param[3])
      
      lambda1 = 0
      
      inf_gene_ind = Cor_ind[1:N1]
      inf_gene = my_genes[inf_gene_ind]
      
      ## Find Diffused genes
      InitialScores = setNames(rep(0,length(my_genes)), my_genes)
      InitialScores[inf_gene]=1
      
      S = page_rank(MyGraph,directed = FALSE,damping = 0.8,
                    personalized = InitialScores)$vector
      S = S/max(S)
      
      s_ind = order(S, decreasing = TRUE)
      s_ind = s_ind[1:N2]
      good_genes = my_genes[s_ind]
      
      X_tr_reduced = X_tr[,good_genes]
      X_val_reduced = X_val[,good_genes]
      
      ## train model
      n_feat = length(good_genes)
      d = rep(1, n_feat)
      GelNet = gelnet(X_tr_reduced, ic50_tr, l1 = lambda1, l2 = lambda2, d = d,
                      P = diag(d), m = rep(0,n_feat), max.iter = 10, eps = 1e-05)
      
      Beta = GelNet[["w"]]
      Beta0 = GelNet[["b"]]
      y_hat = (X_val_reduced %*%  Beta) + Beta0
      y = ic50_val
      
      mse_validation = mean((y - y_hat)^2)
      return(mse_validation)
    }
    
    No_PSOParam = 3                   # N1, N2, and Lambda2
    
    PSO = psoptim(rep(NA,No_PSOParam), fn = FitnessFunc,
                  lower = c(50,1,log(.1)), upper = c(100,70,log(20)),
                  control = list(maxit = 2, s = 2))
    
    ### Training on train + validation
    PSOParam = PSO[["par"]]
    
    N1 = floor(PSOParam[1])
    N2 = floor(PSOParam[2])+N1
    lambda2 = exp(PSOParam[3])
    lambda1 = 0
    
    inf_gene_ind = Cor_ind[1:N1]
    inf_gene = my_genes[inf_gene_ind]
    
    ## Find Diffused genes
    InitialScores = setNames(rep(0,length(my_genes)), my_genes)
    InitialScores[inf_gene]=1
    
    S = page_rank(MyGraph,directed = FALSE,damping = 0.8,
                  personalized = InitialScores)$vector
    S = S/max(S)
    
    s_ind = order(S, decreasing = TRUE)
    s_ind = s_ind[1:N2]
    good_genes = my_genes[s_ind]
    
    X_tr_reduced = X_tr[,good_genes]
    X_val_reduced = X_val[,good_genes]
    X_te_reduced = X_te[,good_genes]
    
    ### Training on train + validation
    n_feat = length(good_genes)
    d = rep(1, n_feat)
    GelNet = gelnet(rbind(X_tr_reduced,X_val_reduced), c(ic50_tr,ic50_val), l1 = lambda1, l2 = lambda2, d = d,
                    P = diag(d), m = rep(0,n_feat), max.iter = 100, eps = 1e-05)
    
    ## Test the model
    Beta = GelNet[["w"]]
    Beta0 = GelNet[["b"]]
    y_hat = (X_te_reduced %*%  Beta) + Beta0
    y = ic50_te
    
    mse_test[j] = mean((y - y_hat)^2)
    corr_test[j] = cor(y,y_hat, method = "pearson")
    
    print(mse_test[j])
    print(corr_test[j])
    
    tempModel = data.frame(N1 = N1, N2 = N2, lambda2 = lambda2, W0 = Beta0)
    Genes = data.frame(rep(0,length(my_genes)))
    row.names(Genes) = my_genes
    
    tempModel = cbind(tempModel,t(Genes))
    tempModel[1,names(Beta)] = Beta
    
    ModelParam = rbind(ModelParam,tempModel)
  }
  
  Test_Result = data.frame(MSE = mse_test, Cor = corr_test)
  print(apply(Test_Result,2,mean))
  print(apply(Test_Result,2,sd))
  
  ## Save File Path
  if (drugi<10) {
    DrugNo = paste0("00",as.character(drugi))
  }else if (drugi<100){
    DrugNo = paste0("0",as.character(drugi))
  }else{
    DrugNo = as.character(drugi)
  }
  FilePath = paste0("Results_InfoGenes/InfoGenes_",DrugNo,".rds")
  
  #save(Test_Result, ModelParam, file = FilePath)
}

