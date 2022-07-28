rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Single_Drug/Infogenes/")

library(STRINGdb)
library(igraph)

##### Load Data
GE = readRDS("Raw_data/expresion_matrix.rds")
sen = readRDS("Raw_data/sensitivity_matrix.rds")
string_db <- STRINGdb$new( version="11.5", species=9606, score_threshold=200, 
                           input_directory="", protocol="http")

full.graph <- string_db$get_graph()
