#                    Created on Wed Aug 11 18:55 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function receives X data to reduce the dimension of 
# genes by finding the transcription factors.

library(decoupleR)
library(doMC)
library(AUCell)
library(ranger)

decoupleR = function(X, method){
  dorothea = get_dorothea(organism = "human", levels = c("A","B","C"))
  X = t(X)
  if (method == "aucell"){
    X_tf = decoupleR::run_aucell(
      mat = X,
      network = dorothea,
      minsize=5,
    )
    }
  
  else if (method == "fgsea"){
    X_tf = run_fgsea(
      mat = X,
      network = dorothea,
      times = 100,
      nproc = 4,
      minsize = 5
    )
    }
  
  else if (method == "gsva"){
    X_tf = run_gsva(
      mat = X,
      network = dorothea,
      method = "gsva",
      minsize = 5
    )
    }
  
  else if (method == "mdt"){
    X_tf = run_mdt(
      mat = X,
      network = dorothea,
      .mor = .data$mor,
      .likelihood = .data$likelihood,
      sparse = FALSE,
      center = FALSE,
      na.rm = FALSE,
      trees = 10,
      min_n = 20,
      nproc = 4,
      minsize = 5
    )}
  
  else if (method == "mlm"){
    X_tf = run_mlm(
      mat = X,
      network = dorothea,
      .mor = .data$mor,
      .likelihood = .data$likelihood,
      sparse = FALSE,
      center = FALSE,
      na.rm = FALSE,
      minsize = 5
    )}
  
  else if (method == "ora"){
    X_tf = run_ora(
      mat = X,
      network = dorothea,
      n_up = ceiling(0.05 * nrow(mat)),
      n_bottom = 0,
      n_background = 20000,
      minsize = 5
    )
    }
  
  else if (method == "udt"){
    X_tf = run_udt(
      mat = X,
      network = dorothea,
      .mor = .data$mor,
      .likelihood = .data$likelihood,
      sparse = FALSE,
      center = FALSE,
      na.rm = FALSE,
      min_n = 20,
      minsize = 5
    )}
  
  else if (method == "ulm"){
    X_tf = run_ulm(
      mat = X,
      network = dorothea,
      .mor = .data$mor,
      .likelihood = .data$likelihood,
      sparse = FALSE,
      center = FALSE,
      na.rm = FALSE,
      minsize = 5
    )}  
  
  else if (method == "viper"){
    X_tf = run_viper(
      mat = X,
      network = dorothea,
      .mor = .data$mor,
      .likelihood = .data$likelihood,
      verbose = FALSE,
      minsize = 5,
      pleiotropy = TRUE,
      eset.filter = FALSE
    )}
  
  else if (method == "wmean"){
    X_tf = run_wmean(
      mat = X,
      network = dorothea,
      .mor = .data$mor,
      .likelihood = .data$likelihood,
      times = 100,
      sparse = TRUE,
      randomize_type = "rows",
      minsize = 5
    )
    }
  
  else if (method == "wsum"){
    X_tf = run_wsum(
      mat = X,
      network = dorothea,
      .mor = .data$mor,
      .likelihood = .data$likelihood,
      times = 100,
      sparse = TRUE,
      randomize_type = "rows",
      minsize = 5 
    )}
  
  else if (method == "consensus"){
    results = decouple(
      mat = X,
      network = dorothea,
      .source = "source",
      .target = "target",
      statistics = c("wmean", "ulm"),
      args = list(
        wmean = list(.mor = "mor", .likelihood = "likelihood"),
        ulm = list(.mor = "mor", .likelihood = "likelihood")
      ),
      consensus_score = FALSE,
      minsize = 0
    )
    X_tf = run_consensus(results)
  }
  
  
  X_tf = data.frame(X_tf)
  X_TF = matrix(0,length(unique(X_tf$condition)), length(unique(X_tf$source)))
  rownames(X_TF) = unique(X_tf$condition) 
  colnames(X_TF) = unique(X_tf$source) 
  
  for (i in rownames(X_TF)){
    for (j in colnames(X_TF)){
      
      if(method == "consensus"){
      X_TF[i,j] = X_tf[X_tf$condition == i & X_tf$source == j,5]
      
      }else{
      X_TF[i,j] = X_tf[X_tf$condition == i & X_tf$source == j,4]
      }
    }
  }
  return(X_TF)
}

