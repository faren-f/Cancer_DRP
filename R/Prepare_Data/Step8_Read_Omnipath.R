rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
library(OmnipathR)
library(tidyr)
library(dnet)
library(gprofiler2)

interactions = import_omnipath_interactions(resources=c("SignaLink3","PhosphoSite",
                                                        "SIGNOR"))

print_interactions(head(interactions))
interactions = import_pathwayextra_interactions(resources=c("BioGRID","STRING"),
                                                organism = 10090)
