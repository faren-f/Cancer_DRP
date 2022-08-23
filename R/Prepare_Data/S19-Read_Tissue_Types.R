rm(list = ls())
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")


# Read Data ---------------------------------------------------------------
cellline_info = read.csv("Raw_data/PRISM/Secondary/secondary-screen-cell-line-info.csv")
cellline_info = cellline_info[complete.cases(cellline_info), ]

sample_ids = unique(cellline_info$depmap_id)
tissue_types = unique(cellline_info$primary_tissue)

sample_tissue_types = matrix(0, length(sample_ids),length(tissue_types))
rownames(sample_tissue_types) = unique(cellline_info$depmap_id)
colnames(sample_tissue_types) = unique(cellline_info$primary_tissue)

for (i in sample_ids){
  for(j in tissue_types)
  
  if(sum(cellline_info$depmap_id==i & cellline_info$primary_tissue==j))
    sample_tissue_types[i,j]=1
  
    else
      sample_tissue_types[i,j]=0
}

GE = readRDS("Processed_Data/S1/expresion_matrix.rds")

sample_tissue_types = sample_tissue_types[intersect(rownames(GE),rownames(sample_tissue_types)),]

saveRDS(sample_tissue_types,"Processed_data/S19/sample_tissue_types.rds")

S = apply(sample_tissue_types,2,sum)
sum(S[S>14])
I = which(S<14)
Other = apply(sample_tissue_types[,I],1,sum)
sample_tissue_types_other = cbind(sample_tissue_types[,-I],Other)

saveRDS(sample_tissue_types_other,"Processed_data/S19/sample_tissue_types_other.rds")


