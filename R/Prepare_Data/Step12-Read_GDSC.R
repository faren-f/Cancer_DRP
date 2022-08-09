rm(list=ls())
library(readr)
library(oligo)
library(affy)
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/Raw_data/GDSC_Raw_data/Cellline_info/GE/")

read.affybatch {affy}
read.affybatch("data.csv", filenames = character(0),
               phenoData = new("AnnotatedDataFrame"),
               compress = getOption("BioC")$affy$compress.cel)
               
ReadAffy(filenames=character(0),
         widget=getOption("BioC")$affy$use.widgets,
         compress=getOption("BioC")$affy$compress.cel,
         celfile.path="E-MTAB-3610.raw.1/")
