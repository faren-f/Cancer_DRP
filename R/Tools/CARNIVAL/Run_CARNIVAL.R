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


#-----------------------------------------------------------
load(file = system.file("toy_perturbations_ex1.RData",
                        package="CARNIVAL"))
load(file = system.file("toy_measurements_ex1.RData",
                        package="CARNIVAL"))
load(file = system.file("toy_network_ex1.RData",
                        package="CARNIVAL"))

res1 = runCARNIVAL(inputObj = toy_perturbations_ex1,
                   measObj = toy_measurements_ex1,
                   netObj = toy_network_ex1,
                   solverPath = "/Users/faren/Softwares/CPLEX_Studio_Community221/cplex/bin/x86-64_osx/cplex",
                   solver = 'cplex')

res1$weightedSIF ##see @return
res1$nodesAttributes ## see @return
res1$sifAll ## see @return
res1$attributesAll ## see @return


