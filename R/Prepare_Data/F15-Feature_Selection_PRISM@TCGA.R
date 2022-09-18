#                  Created on Wed Aug 20 15:44 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: 


source("F9-decoupleR.R")
source("F4-DoRothEA.R")
source("F5-Progeny.R")
 
Feature_Selection = function(selected_features,GE,GE_test){
  if (prod(selected_features == "")){
    writeLines("selected_features is empty!\nEnter your desired features")
    omics = c()
    index = c()
    Omics = list(omics,index)
    
  }else if (prod(selected_features == "Whole_genes")){
    omics = GE
    omics_test = GE_test
    
    index = rep(1,ncol(omics))
    Omics = list(omics,index,omics_test)
    
  }else if (prod(selected_features == "Landmark_genes")){
    l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
    I_G = intersect(l1000_genes,colnames(GE_test))
    
    omics = GE[,I_G]
    omics_test = GE_test[,I_G]
    
    index = rep(1,ncol(omics))
    Omics = list(omics,index,omics_test)
    
  }else if (prod(selected_features == "TF_DoRothEA")){
    omics = DoRothEA(GE)
    omics_test = DoRothEA(GE_test)
    
    index = rep(1,ncol(omics))
    Omics = list(omics,index,omics_test)
    
  }else if (prod(selected_features == "Tissue_types")){
    omics = readRDS("Processed_data/S19/sample_tissue_types.rds")
    index = rep(1,ncol(omics))
    Omics = list(omics,index)    
    
  }else if (prod(selected_features == "TF_decoupleR")){
    omics = decoupleR(X = GE, method = "gsva")
    omics_test = decoupleR(X = GE_test, method = "gsva")
    
    index = rep(1,ncol(omics))
    Omics = list(omics,index,omics_test)
    
  }else if (prod(selected_features == "Pathway_genes")){
    omics = Pathway_genes(X = GE, drug=)
    index = rep(1,ncol(omics))
    Omics = list(omics,index)
    
  }else if (prod(selected_features == "progeny")){
    omics = Progeny_pw_act(GE)
    omics_test = Progeny_pw_act(GE_test)
    
    index = rep(1,ncol(omics))
    Omics = list(omics,index,omics_test)
    
  }else if (prod(selected_features == c("Landmark_genes","TF_DoRothEA"))){
    l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
    I_G = intersect(l1000_genes,colnames(GE_test))
    
    O1 = GE[,I_G]
    O1_test = GE_test[,I_G]
    
    O2 = DoRothEA(GE)
    O2_test = DoRothEA(GE_test) 
    
    omics = cbind(O1,O2)
    omics_test = cbind(O1_test,O2_test)
    
    colnames(omics) = 1:ncol(omics)
    colnames(omics_test) = 1:ncol(omics_test)
    
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)))
    Omics = list(omics,index,omics_test)
    
  }else if (prod(selected_features == c("Landmark_genes","Tissue_types"))){
    l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
    O1 = GE[,colnames(GE)%in%l1000_genes]
    O2 = readRDS("Processed_data/S19/sample_tissue_types.rds")
    omics = cbind(O1,O2)
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)))
    Omics = list(omics,index)
    
  }else if (prod(selected_features == c("Landmark_genes","TF_decoupleR"))){
    l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
    I_G = intersect(l1000_genes,colnames(GE_test))
    
    O1 = GE[,I_G]
    O1_test = GE_test[,I_G]
    
    O2 = decoupleR(X = GE, method = "gsva")
    O2_test = decoupleR(X = GE_test, method = "gsva")

    
    omics = cbind(O1,O2)
    omics_test = cbind(O1_test,O2_test)
    
    colnames(omics) = 1:ncol(omics)
    colnames(omics_test) = 1:ncol(omics_test)
    
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)))
    Omics = list(omics,index,omics_test)
    
  }else if (prod(selected_features == c("TF_DoRothEA","Tissue_types"))){
    
    O1 = DoRothEA(GE)
    O1_test = DoRothEA(GE_test)
    
    ##incompelete
    O2 = readRDS("Processed_data/S19/sample_tissue_types.rds")
    
    
    omics = cbind(O1,O2)
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)))
    Omics = list(omics,index)
    
    
  }else if (prod(selected_features == c("Landmark_genes","progeny"))){
    l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
    I_G = intersect(l1000_genes,colnames(GE_test))
    
    O1 = GE[,I_G]
    O1_test = GE_test[,I_G]
    
    O2 = Progeny_pw_act(GE)
    O2_test = Progeny_pw_act(GE_test)
    
    omics = cbind(O1,O2)
    omics_test = cbind(O1_test,O2_test)
    
    colnames(omics) = 1:ncol(omics)
    colnames(omics_test) = 1:ncol(omics_test)
    
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)))
    Omics = list(omics,index,omics_test)
    
  }else if (prod(selected_features == c("TF_DoRothEA","progeny"))){
    O1 = DoRothEA(GE)
    O1_test = DoRothEA(GE_test)
    
    O2 = Progeny_pw_act(GE)
    O2_test = Progeny_pw_act(GE_test)
    
    omics = cbind(O1,O2)
    omics_test = cbind(O1_test,O2_test)
    
    colnames(omics) = 1:ncol(omics)
    colnames(omics_test) = 1:ncol(omics_test)
    
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)))
    Omics = list(omics,index,omics_test)
    
  }else if (prod(selected_features == c("TF_decoupleR","progeny"))){
    
    
    O1 = decoupleR(X = GE, method = "gsva")
    O1_test = decoupleR(X = GE_test, method = "gsva")
    
    O2 = Progeny_pw_act(GE)
    O2_test = Progeny_pw_act(GE_test)
    
    omics = cbind(O1,O2)
    omics_test = cbind(O1_test,O2_test)
    
    colnames(omics) = 1:ncol(omics)
    colnames(omics_test) = 1:ncol(omics_test)
    
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)))
    Omics = list(omics,index,omics_test)
    
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
    I_G = intersect(l1000_genes,colnames(GE_test))
    
    O1 = GE[,I_G]
    O1_test = GE_test[,I_G]
    
    O2 = DoRothEA(GE)
    O2_test = DoRothEA(GE_test) 
    
    O3 = decoupleR(X = GE, method = "gsva")
    O3_test = decoupleR(X = GE_test, method = "gsva")
    
    omics = cbind(O1,O2,O3)
    omics_test = cbind(O1_test,O2_test,O3_test)
    
    colnames(omics) = 1:ncol(omics)
    colnames(omics_test) = 1:ncol(omics_test)
    
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)),rep(3,ncol(O3)))
    Omics = list(omics,index,omics_test)
    
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


