#                  Created on Wed Aug 20 15:44 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: 

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
conv_table = readRDS("Processed_data/S7/biomart_conversion_table.rds")

GE = readRDS("Processed_Data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")


### for cisplatin
# GE = readRDS("Processed_Data/S1/expresion_matrix.rds")
# cisplatin_genes = readRDS("Processed_data/S26/cisplatin_gene_pathways.rds")
# I = intersect(colnames(GE),cisplatin_genes$hgnc_symbol)
#GE = GE[,a]
##########


source("F9-decoupleR.R")
source("F5-Progeny.R")

Feature_Selection = function(selected_features){
  Omics = list()
  if (prod(selected_features == "")){
    writeLines("selected_features is empty!\nEnter your desired features")
    omics = c()
    index = c()
    Omics = list(omics,index)
  
  }else if (prod(selected_features == "Whole_genes")){
    omics = GE
    index = rep(1,ncol(omics))
    Omics = list(omics,index)
    
  }else if (prod(selected_features == "Landmark_genes")){
    l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
    omics = GE[,colnames(GE)%in%l1000_genes]
    index = rep(1,ncol(omics))
    Omics = list(omics,index)
    
  }else if (prod(selected_features == "TF_DoRothEA")){
    omics = readRDS("Processed_data/S14/DoRothEA_TF.rds")
    index = rep(1,ncol(omics))
    Omics = list(omics,index)
    
  }else if (prod(selected_features == "Tissue_types")){
    omics = readRDS("Processed_data/S19/sample_tissue_types.rds")
    index = rep(1,ncol(omics))
    Omics = list(omics,index)    
    
  }else if (prod(selected_features == "TF_decoupleR")){
    omics = decoupleR(X = GE, method = "gsva")
    index = rep(1,ncol(omics))
    Omics = list(omics,index)
  
  }else if (prod(selected_features == "progeny")){
    omics = Progeny_pw_act(X = GE)
    index = rep(1,ncol(omics))
    Omics = list(omics,index)
  }else if (prod(selected_features == c("Landmark_genes","TF_DoRothEA"))){
    l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
    O1 = GE[,colnames(GE)%in%l1000_genes]
    O2 = readRDS("Processed_data/S14/DoRothEA_TF.rds")

    omics = cbind(O1,O2)
    colnames(omics) = 1:ncol(omics)
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)))
    Omics = list(omics,index)
    
  }else if (prod(selected_features == c("Landmark_genes","Tissue_types"))){
    l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
    O1 = GE[,colnames(GE)%in%l1000_genes]
    O2 = readRDS("Processed_data/S19/sample_tissue_types.rds")
    omics = cbind(O1,O2)
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)))
    Omics = list(omics,index)
    
  }else if (prod(selected_features == c("Landmark_genes","TF_decoupleR"))){
    l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
    O1 = GE[,colnames(GE)%in%l1000_genes]
    O2 = decoupleR(X = GE, method = "gsva")

    omics = cbind(O1,O2)
    colnames(omics) = 1:ncol(omics)
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)))
    Omics = list(omics,index)
    
  }else if (prod(selected_features == c("TF_DoRothEA","Tissue_types"))){
    O1 = readRDS("Processed_data/S14/DoRothEA_TF.rds")
    O2 = readRDS("Processed_data/S19/sample_tissue_types.rds")
    omics = cbind(O1,O2)
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)))
    Omics = list(omics,index)
    
  }else if (prod(selected_features == c("Landmark_genes","TF_DoRothEA","Tissue_types"))){
    l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
    O1 = GE[,colnames(GE)%in%l1000_genes]
    O2 = readRDS("Processed_data/S14/DoRothEA_TF.rds")
    O3 = readRDS("Processed_data/S19/sample_tissue_types.rds")
    
    omics = cbind(O1,O2,O3)
    colnames(omics) = 1:ncol(omics)
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)),rep(3,ncol(O3)))
    Omics = list(omics,index)
    
  }else if (prod(selected_features == c("Landmark_genes","TF_DoRothEA","TF_decoupleR"))){
    l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
    O1 = GE[,colnames(GE)%in%l1000_genes]
    O2 = readRDS("Processed_data/S14/DoRothEA_TF.rds")
    O3 = decoupleR(X = GE, method = "gsva")

    omics = cbind(O1,O2,O3)
    colnames(omics) = 1:ncol(omics)
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)),rep(3,ncol(O3)))
    Omics = list(omics,index)    
    
  }else if (prod(selected_features == c("Landmark_genes","Tissue_types","TF_decoupleR"))){
    l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
    O1 = GE[,colnames(GE)%in%l1000_genes]
    O2 = readRDS("Processed_data/S19/sample_tissue_types.rds")
    O3 = decoupleR(X = GE, method = "gsva")

    omics = cbind(O1,O2,O3)
    colnames(omics) = 1:ncol(omics)
    
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)),rep(3,ncol(O3)))
    Omics = list(omics,index)    
    
  }else if (prod(selected_features == c("Landmark_genes","Tissue_types","TF_DoRothEA","TF_decoupleR"))){
    
    l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
    O1 = GE[,colnames(GE)%in%l1000_genes]
    O2 = readRDS("Processed_data/S19/sample_tissue_types.rds")
    O3 = readRDS("Processed_data/S14/DoRothEA_TF.rds")
    O4 = decoupleR(X = GE, method = "gsva")
    
    omics = cbind(O1,O2,O3,O4)
    colnames(omics) = 1:ncol(omics)
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)),rep(3,ncol(O3)),rep(4,ncol(O4)))
    Omics = list(omics,index)    
  }
return(Omics)
}


