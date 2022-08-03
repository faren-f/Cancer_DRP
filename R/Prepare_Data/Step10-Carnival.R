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

setwd("~/Desktop/Cancer_DRP/R/Single_Drug/RF/")

GE = readRDS("Raw_data/expresion_matrix.rds")
sen = readRDS("Raw_data/sensitivity_matrix.rds")



#-----------------------------------------------------------
load(file = system.file("toy_perturbations_ex1.RData",
                        package="CARNIVAL"))
load(file = system.file("toy_measurements_ex1.RData",
                        package="CARNIVAL"))
load(file = system.file("toy_network_ex1.RData",
                        package="CARNIVAL"))

## lpSolve
#res1 = runCARNIVAL(inputObj = toy_perturbations_ex1,
#                    measObj = toy_measurements_ex1,
#                    netObj = toy_network_ex1,
#                    solver = 'lpSolve')

#res1$weightedSIF ##see @return
#res1$nodesAttributes ## see @return
#res1$sifAll ## see @return
#res1$attributesAll ## see @return

## Examples for cbc and cplex are commented out because these solvers are not part of R environment
## and need to be installed separately
##
## cbc
## res2 = runCARNIVAL(inputObj = toy_perturbations_ex1,
##                    measObj = toy_measurements_ex1,
##                    netObj = toy_network_ex1,
##                    solver = 'cbc')
##
## res2$weightedSIF ##see @return
## res2$nodesAttributes ## see @return
## res2$sifAll ## see @return
## res2$attributesAll ## see @return
##
## cplex
## res3 = runCARNIVAL(inputObj = toy_perturbations_ex1,
##                    measObj = toy_measurements_ex1,
##                    netObj = toy_network_ex1,
##                    solver = 'cplex')
##
## res3$weightedSIF ##see @return
## res3$nodesAttributes ## see @return
## res3$sifAll ## see @return
## res3$attributesAll ## see @return


runCARNIVAL(
  inputObj = NULL,
  measObj = measObj,
  netObj = netObj,
  weightObj = NULL,
  solverPath = NULL,
  solver = c("lpSolve", "cplex", "cbc", "gurobi"),
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





runVanillaCarnival(
  perturbations,
  measurements,
  priorKnowledgeNetwork,
  weights = NULL,
  carnivalOptions = defaultLpSolveCarnivalOptions()
)

checkCarnivalOptions(carnivalOptions)





