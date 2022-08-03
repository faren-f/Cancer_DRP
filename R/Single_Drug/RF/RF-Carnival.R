rm(list=ls())

library(CARNIVAL)
library(readr)
library(dplyr)
library(lpSolve)
library(stringr)
library(igraph)
library(tibble)
library(tidyr)
library(rjson)
library(rmarkdown)
library(RefManageR)
library(BiocStyle)
library(covr)
library(knitr)
library(testthat)
library(sessioninfo)


library(dorothea)
library(OmnipathR)


setwd("~/Desktop/Cancer_DRP/R/Single_Drug/RF/")

GE = readRDS("Raw_data/expresion_matrix.rds")
sen = readRDS("Raw_data/sensitivity_matrix.rds")


i = 325                            # drug number
GE = GE[!is.na(sen[,i]),]
sen_i = sen[!is.na(sen[,i]),i]

# Finding differentially expressed genes (DEGs)
Cor_GE_sen = apply(GE,2,function(x){return(abs(cor(x,sen_i)))})
idx_Cor = order(Cor_GE_sen, decreasing = TRUE)

GE = GE[,1:200]
GE = t(GE)              # input: rows are genes and columns are cell lines 

# Finding transcription activities using dorothea
data(dorothea_hs, package = "dorothea")

tf_activities <- run_viper(GE, dorothea_hs,
                           options =  list(method = "scale", minsize = 4,
                                           eset.filter = FALSE, cores = 1,
                                           verbose = FALSE))
TFs = rownames(tf_activities)
measurements = rep(1, length(TFs))
names(measurements) = TFs


# Download protein-protein interactions
interactions = import_omnipath_interactions() %>% as_tibble()
net = interactions[,c(3,6,4)]

colnames(net) = c("source","interaction","target")

n = replace(net$interaction, net$interaction==0,-1)
n = ifelse(net$interaction==0,-1,1)
net[,2]=n

priorKnowledgeNetwork = net

drug_targets = c("ABL1","SRC","EPHA2","YES1","KITLG")
perturbations = rep(1,length(drug_targets))
names(perturbations) = drug_targets



runCARNIVAL(
  inputObj = perturbations,
  measObj = measurements,
  netObj = priorKnowledgeNetwork,
  weightObj = NULL,
  solverPath = "/Users/faren/",
  solver = "cplex",
  timelimit = 3600,
  mipGAP = 0.05,
  poolrelGAP = 1e-04,
  limitPop = 500,
  poolCap = 100,
  poolIntensity = 4,
  poolReplace = 2,
  alphaWeight = 1,
  betaWeight = 0.2,
  threads = 0,
  cleanTmpFiles = TRUE,
  keepLPFiles = TRUE,
  clonelog = -1,
  dir_name = getwd()
)

#https://github.com/saezlab/CARNIVAL/issues/51









runVanillaCarnival(
  perturbations,
  measurements,
  priorKnowledgeNetwork,
  weights = NULL,
  carnivalOptions = defaultLpSolveCarnivalOptions()
)

checkCarnivalOptions(carnivalOptions)




