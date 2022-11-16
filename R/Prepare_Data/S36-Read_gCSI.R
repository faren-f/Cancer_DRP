rm(list=ls())
library(PharmacoGx)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

GDSCv2 = readRDS("Processed_data/S34/gCSI_PharmacoGx.rds")

CellName = phenoInfo(GDSCv2,"rna")                                  
CellInfo = cellInfo(GDSCv2)                                         
CellInfo = CellInfo[CellName$Characteristics.cell.line.,]        

GeneName = featureInfo(GDSCv2, "rna")                               
Expr = molecularProfiles(GDSCv2)                         

rownames(Expr) = CellName$Characteristics.cell.line.
colnames(Expr) = GeneName$Symbol

sen = summarizeSensitivityProfiles(pSet = GDSCv2,
                                   sensitivity.measure="ic50_published", 
                                   summary.stat="median") 



