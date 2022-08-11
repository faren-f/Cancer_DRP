#                    Created on Wed Aug 10 12:04 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function receives X data to reduce the dimension of 
# genes by finding the transcription factors.
rm(list=ls())

library(decoupleR)
library(doMC)
library(AUCell)
library(ranger)

show_methods()
all_methods = c("aucell", "fgsea", "gsva", "mdt", "mlm", "ora", "udt", "ulm", 
                "viper", "wmean", "wsum", "consensus")
dorothea = get_dorothea(organism = "human", levels = c("A","B","C"))

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
X = readRDS("Processed_Data/Step1/expresion_matrix.rds")
X = t(X)
method = "wmean"

decoupleR_DR = function(X, dorothea, method){
  
  if (method == "aucell"){
    X_TF = run_aucell(
      mat = X,
      network = dorothea,
      minsize=5,
    )}
  
  else if (method == "fgsea"){
    X_TF = run_fgsea(
      mat = X,
      network = dorothea,
      times = 100,
      nproc = 4,
      minsize = 5
    )}
    
  else if (method == "gsva"){
    X_TF = run_gsva(
      mat = X,
      network = dorothea,
      method = "gsva",
      minsize = 5
    )}
    
  else if (method == "mdt"){
    X_TF = run_mdt(
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
    X_TF = run_mlm(
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
    X_TF = run_ora(
    mat = X,
    network = dorothea,
    n_up = ceiling(0.05 * nrow(mat)),
    n_bottom = 0,
    n_background = 20000,
    minsize = 4
    )}
  
  else if (method == "udt"){
    X_TF = run_udt(
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
    X_TF = run_ulm(
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
    X_TF = run_viper(
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
    X_TF = run_wmean(
      mat = X,
      network = dorothea,
      .mor = .data$mor,
      .likelihood = .data$likelihood,
      times = 100,
      sparse = TRUE,
      randomize_type = "rows",
      minsize = 5
    )}
  
  else if (method == "wsum"){
    X_TF = run_wsum(
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
    X_TF = run_consensus(results)
  }
  return(X_TF)
}



X_TF = decoupleR_DR(X, dorothea, method)

saveRDS(X_TF,"Processed_data/Step13/TF_wmean.rds")

