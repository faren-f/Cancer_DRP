rm(list=ls())
library(PharmacoGx)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
GDSC = readRDS("Processed_data/S38/GDSC_PharmacoGx.rds")

CellName = phenoInfo(GDSC,"rna") 
GeneName = featureInfo(GDSC, "rna")                               
#DrugInfo = drugInfo(GDSC)

#CellInfo = cellInfo(GDSC)                                         
#CellInfo = CellInfo[CellName$Characteristics.cell.line.,]        
#sen = sensitivityInfo(GDSC)

GE = t(molecularProfiles(GDSC, "rna"))                         
rownames(GE) = CellName$Characteristics.cell.line.
colnames(GE) = GeneName$Symbol

sen = summarizeSensitivityProfiles(GDSC, sensitivity.measure='ic50_recomputed')                                  
sen = t(sen)
sen = log(sen)

I = intersect(rownames(GE),rownames(sen))
GE = GE[I,]
sen = sen[I,]

hist(sen,1000)
hist(sen, xlim = c(200,1000),ylim = c(0,2),300)

