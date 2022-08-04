rm(list=ls())

library(CARNIVAL)
library(dplyr)
library(tibble)
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

GE = GE[,1:600]
GE = t(GE)              # input: rows are genes and columns are cell lines 

# Finding transcription activities using dorothea
data(dorothea_hs, package = "dorothea")

tf_activities <- run_viper(GE, dorothea_hs,
                           options =  list(method = "scale", minsize = 4,
                                           eset.filter = FALSE, cores = 1,
                                           verbose = FALSE))
TFs = rownames(tf_activities)



# Download protein-protein interactions
interactions = import_omnipath_interactions() %>% as_tibble()
net = interactions[,c(3,6,4)]

a = mutate(net,in_OP_source = net$source_genesymbol %in% TFs)
a = mutate(a,in_OP_target = a$target_genesymbol %in% TFs)
net = net[which(a$in_OP_source&a$in_OP_target),]
colnames(net) = c("source","interaction","target")

s = unique(c(net$source,net$target))
TFs = intersect(TFs,s)


n = replace(net$interaction, net$interaction==0,-1)
n = ifelse(net$interaction==0,-1,1)
net[,2]=n

#drug_targets = c("ABL1","SRC","EPHA2","YES1","KITLG")
#perturbations = rep(1,length(drug_targets))
#names(perturbations) = drug_targets


measurements = rep(1, length(TFs))
names(measurements) = TFs

priorKnowledgeNetwork = net


res = runCARNIVAL(
  inputObj = NULL,
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
res3attributesAll ## see @return






load(file = system.file("toy_perturbations_ex1.RData",
                        package="CARNIVAL"))
load(file = system.file("toy_measurements_ex1.RData",
                        package="CARNIVAL"))
load(file = system.file("toy_network_ex1.RData",
                        package="CARNIVAL"))

res3 = runCARNIVAL(inputObj = toy_perturbations_ex1,
                   measObj = toy_measurements_ex1,
                   netObj = toy_network_ex1,
                   solverPath = "/Users/faren/Softwares/CPLEX_Studio_Community221/cplex/bin/x86-64_osx/cplex",
                   solver = 'cplex')

res3$weightedSIF ##see @return
res3$nodesAttributes ## see @return
res3$sifAll ## see @return
res3$attributesAll ## see @return












