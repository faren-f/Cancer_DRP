rm(list=ls())
library(PharmacoGx)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
GDSC = readRDS("Processed_data/S38/GDSC_PharmacoGx.rds")

CellName = phenoInfo(GDSC,"rna")                                  
CellInfo = cellInfo(GDSC)                                         
CellInfo = CellInfo[CellName$Characteristics.cell.line.,]        

GeneName = featureInfo(GDSC, "rna")                               
GE = t(molecularProfiles(GDSC, "rna"))                         

rownames(GE) = CellName$Characteristics.cell.line.
colnames(GE) = GeneName$Symbol

Sen = summarizeSensitivityProfiles(GDSC, sensitivity.measure='ic50_recomputed')                                  
Sen = t(Sen)

# sen = data.frame(sensitivityRaw(GDSC))
# sen = treatmentResponse(GDSC)
# s = data.frame(sensitivityMeasures(GDSC))
# DrugInfo = drugInfo(GDSC)
# sen = sensitivityInfo(GDSC)

I = intersect(rownames(GE),rownames(Sen))
GE = GE[I,]
Sen = Sen[I,]







