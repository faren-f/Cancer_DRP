rm(list=ls())

library(CARNIVAL)
library(dplyr)
library(tibble)
library(dorothea)
library(OmnipathR)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

GE = readRDS("Processed_Data/Step1/expresion_matrix.rds")
sen = readRDS("Processed_Data/Step1/sensitivity_matrix.rds")

i = 325                            # drug number
GE = GE[!is.na(sen[,i]),]
sen_i = sen[!is.na(sen[,i]),i]

# Finding differentially expressed genes (DEGs)
Cor_GE_sen = apply(GE,2,function(x){return(abs(cor(x,sen_i)))})
idx_Cor = order(Cor_GE_sen, decreasing = TRUE)

GE = GE[,1:600]
GE = t(GE)              # input: rows are genes and columns are cell lines 

# Finding transcription activities using DoRothEA
data(dorothea_hs, package = "dorothea")
tf_activities <- run_viper(GE, dorothea_hs,
                           options =  list(method = "scale", minsize = 4,
                                           eset.filter = FALSE, cores = 1,
                                           verbose = FALSE))
TFs = rownames(tf_activities)

# Download protein-protein interactions
interactions = import_omnipath_interactions() %>% as_tibble()
net = interactions[,c(3,6,4)]

# Finding overlaps between genes that are available in Omnipath interaction and TFs
net = net[net$source_genesymbol %in% TFs & net$target_genesymbol %in% TFs,]
colnames(net) = c("source","interaction","target")
n = replace(net$interaction, net$interaction==0,-1)
n = ifelse(net$interaction==0,-1,1)
net[,2]=n

TFs = intersect(TFs,unique(c(net$source,net$target)))

# Drug targets
drug_targets = c("ABL1","SRC","EPHA2","YES1","KITLG")
drug_targets = intersect(drug_targets,TFs)

perturbations = rep(1,length(drug_targets))
names(perturbations) = drug_targets

if (length(perturbations)==0){
  perturbations = NULL
}else{
  perturbations = perturbations
}

measurements = rep(1, length(TFs))
names(measurements) = TFs

priorKnowledgeNetwork = net

res = runCARNIVAL(
  inputObj = perturbations,
  measObj = measurements,
  netObj = priorKnowledgeNetwork,
  weightObj = NULL,
  solverPath = "/Users/faren/Softwares/CPLEX_Studio_Community221/cplex/bin/x86-64_osx/cplex",
  solver = "cplex",
  timelimit = 7200,
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

res$weightedSIF ##see @return
res$nodesAttributes ## see @return
res$sifAll ## see @return
res$attributesAll ## see @return

