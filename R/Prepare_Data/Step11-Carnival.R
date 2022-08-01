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


runVanillaCarnival(
  perturbations,
  measurements,
  priorKnowledgeNetwork,
  weights = NULL,
  carnivalOptions = defaultLpSolveCarnivalOptions()
)

checkCarnivalOptions(carnivalOptions)


https://rdrr.io/github/saezlab/CARNIVAL/api/ 



