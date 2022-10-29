rm(list=ls())
library(PharmacoGx)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

GDSCv2 = readRDS("Processed_data/S38/GDSC_v2_PharmacoGx.rds")

CellName = phenoInfo(GDSCv2,"rna")                                  
CellInfo = cellInfo(GDSCv2)                                         
CellInfo = CellInfo[CellName$Characteristics.cell.line.,]        

GeneName = featureInfo(GDSCv2, "rna")                               
Expr = t(molecularProfiles(GDSCv2, "rna"))                         

rownames(Expr) = CellName$Characteristics.cell.line.
colnames(Expr) = GeneName$Symbol


sensitivity = summarizeSensitivityProfiles(GDSCv2, sensitivity.measure='ic50_recomputed')                                  


sen = data.frame(sensitivityRaw(GDSCv2))
s = data.frame(sensitivityMeasures(GDSCv2))
a = treatmentResponse(GDSCv2)

DrugInfo = drugInfo(GDSCv2)
sen = sensitivityInfo(GDSCv2)



