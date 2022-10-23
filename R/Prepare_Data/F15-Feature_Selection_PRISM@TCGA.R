#                  Created on Wed Aug 20 15:44 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: 


source("F9-decoupleR.R")
source("F4-DoRothEA.R")
source("F5-Progeny.R")
 
Feature_Selection_PRISM_TCGA = function(selected_features,Xtrain,Xtest){
  if (prod(selected_features == "")){
    writeLines("selected_features is empty!\nEnter your desired features")
    omics_train = c()
    omics_test = c()
    index = c()
    Omics = list(omics_train,index,omics_test)
    
  }else if (prod(selected_features == "Whole_genes")){
    omics_train = Xtrain
    omics_test = Xtest
    
    index = rep(1,ncol(omics_train))
    Omics = list(omics_train,index,omics_test)
    
  }else if (prod(selected_features == "Landmark_genes")){
    l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
    I_G = intersect(l1000_genes,colnames(Xtest))
    
    omics_train = Xtrain[,I_G]
    omics_test = Xtest[,I_G]
    
    index = rep(1,ncol(omics_train))
    Omics = list(omics_train,index,omics_test)
    
  }else if (prod(selected_features == "TF_DoRothEA")){
    omics_train = DoRothEA(Xtrain)
    omics_test = DoRothEA(Xtest)
    
    index = rep(1,ncol(omics_train))
    Omics = list(omics_train,index,omics_test)
    
  }else if (prod(selected_features == "Tissue_types")){
    omics_train = readRDS("Processed_data/S19/sample_tissue_types.rds")
    index = rep(1,ncol(omics_train))
    Omics = list(omics_train,index)    
    
  }else if (prod(selected_features == "TF_decoupleR")){
    omics_train = decoupleR(X = Xtrain, method = "ulm")
    omics_test = decoupleR(X = Xtest, method = "ulm")
    
    index = rep(1,ncol(omics_train))
    Omics = list(omics_train,index,omics_test)
    
  }else if (prod(selected_features == "Pathway_genes")){
    omics_train = Pathway_genes(X = Xtrain, drug=)
    index = rep(1,ncol(omics_train))
    Omics = list(omics_train,index)
    
  }else if (prod(selected_features == "progeny")){
    omics_train = Progeny_pw_act(Xtrain)
    omics_test = Progeny_pw_act(Xtest)
    
    index = rep(1,ncol(omics_train))
    Omics = list(omics_train,index,omics_test)
    
  }else if (prod(selected_features == "OncoKB")){
    OncoKB_genes = readRDS("Processed_Data/S32/OncoKB_genes.rds")
    I_G = intersect(OncoKB_genes,colnames(Xtest))
    
    omics_train = Xtrain[,I_G]
    omics_test = Xtest[,I_G]
    
    index = rep(1,ncol(omics_train))
    Omics = list(omics_train,index,omics_test)
    
  }else if (prod(selected_features == "CancerGenes")){
    CancerGenesGDSC = readRDS("Processed_data/S33/CancerGenesGDSC.rds")
    I_G = intersect(CancerGenesGDSC,colnames(Xtest))
    
    omics_train = Xtrain[,I_G]
    omics_test = Xtest[,I_G]
    
    index = rep(1,ncol(omics_train))
    Omics = list(omics_train,index,omics_test)
    
    
  }else if (prod(selected_features == c("Landmark_OncoKB_genes"))){
    l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
    OncoKB = readRDS("Processed_Data/S32/OncoKB_genes.rds")
    Landmark_OncoKB_genes = intersect(OncoKB, l1000_genes)

    I_G = intersect(Landmark_OncoKB_genes,colnames(Xtest))
    
    omics_train = Xtrain[,I_G]
    omics_test = Xtest[,I_G]
    
    index = rep(1,ncol(omics_train))
    Omics = list(omics_train,index,omics_test)
    
  }else if (prod(selected_features == c("Landmark_genes","OncoKB"))){
    l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
    OncoKB = readRDS("Processed_Data/S32/OncoKB_genes.rds")
    Landmark_OncoKB_genes = c(OncoKB, l1000_genes)
    Landmark_OncoKB_genes = Landmark_OncoKB_genes[!duplicated(Landmark_OncoKB_genes)] 
    I_G = intersect(Landmark_OncoKB_genes,colnames(Xtest))
    
    omics_train = Xtrain[,I_G]
    omics_test = Xtest[,I_G]
    
    index = rep(1,ncol(omics_train))
    Omics = list(omics_train,index,omics_test)
        
  }else if (prod(selected_features == "OncoKB_oncogenes")){
    OncoKB_oncogenes = readRDS("Processed_Data/S32/OncoKB_oncogenes.rds")
    I_G = intersect(OncoKB_oncogenes,colnames(Xtest))
    
    omics_train = Xtrain[,I_G]
    omics_test = Xtest[,I_G]
    
    index = rep(1,ncol(omics_train))
    Omics = list(omics_train,index,omics_test)
    
    
  }else if (prod(selected_features == c("Landmark_genes","TF_DoRothEA"))){
    l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
    I_G = intersect(l1000_genes,colnames(Xtest))
    
    O1 = Xtrain[,I_G]
    O1_test = Xtest[,I_G]
    
    O2 = DoRothEA(Xtrain)
    O2_test = DoRothEA(Xtest) 
    
    omics_train = cbind(O1,O2)
    omics_test = cbind(O1_test,O2_test)
    
    colnames(omics_train) = 1:ncol(omics_train)
    colnames(omics_test) = 1:ncol(omics_test)
    
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)))
    Omics = list(omics_train,index,omics_test)
    
  }else if (prod(selected_features == c("Landmark_genes","Tissue_types"))){
    l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
    O1 = Xtrain[,colnames(Xtrain)%in%l1000_genes]
    O2 = readRDS("Processed_data/S19/sample_tissue_types.rds")
    omics_train = cbind(O1,O2)
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)))
    Omics = list(omics_train,index)
    
  }else if (prod(selected_features == c("Landmark_genes","TF_decoupleR"))){
    l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
    I_G = intersect(l1000_genes,colnames(Xtest))
    
    O1 = Xtrain[,I_G]
    O1_test = Xtest[,I_G]
    
    O2 = decoupleR(X = Xtrain, method = "gsva")
    O2_test = decoupleR(X = Xtest, method = "gsva")

    
    omics_train = cbind(O1,O2)
    omics_test = cbind(O1_test,O2_test)
    
    colnames(omics_train) = 1:ncol(omics_train)
    colnames(omics_test) = 1:ncol(omics_test)
    
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)))
    Omics = list(omics_train,index,omics_test)
    
  }else if (prod(selected_features == c("TF_DoRothEA","Tissue_types"))){
    
    O1 = DoRothEA(Xtrain)
    O1_test = DoRothEA(Xtest)
    
    ##incompelete
    O2 = readRDS("Processed_data/S19/sample_tissue_types.rds")
    
    
    omics_train = cbind(O1,O2)
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)))
    Omics = list(omics_train,index)
    
    
  }else if (prod(selected_features == c("Landmark_genes","progeny"))){
    l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
    I_G = intersect(l1000_genes,colnames(Xtest))
    
    O1 = Xtrain[,I_G]
    O1_test = Xtest[,I_G]
    
    O2 = Progeny_pw_act(Xtrain)
    O2_test = Progeny_pw_act(Xtest)
    
    omics_train = cbind(O1,O2)
    omics_test = cbind(O1_test,O2_test)
    
    colnames(omics_train) = 1:ncol(omics_train)
    colnames(omics_test) = 1:ncol(omics_test)
    
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)))
    Omics = list(omics_train,index,omics_test)
    
  }else if (prod(selected_features == c("TF_DoRothEA","progeny"))){
    O1 = DoRothEA(Xtrain)
    O1_test = DoRothEA(Xtest)
    
    O2 = Progeny_pw_act(Xtrain)
    O2_test = Progeny_pw_act(Xtest)
    
    omics_train = cbind(O1,O2)
    omics_test = cbind(O1_test,O2_test)
    
    colnames(omics_train) = 1:ncol(omics_train)
    colnames(omics_test) = 1:ncol(omics_test)
    
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)))
    Omics = list(omics_train,index,omics_test)
    
  }else if (prod(selected_features == c("TF_decoupleR","progeny"))){
    
    
    O1 = decoupleR(X = Xtrain, method = "gsva")
    O1_test = decoupleR(X = Xtest, method = "gsva")
    
    O2 = Progeny_pw_act(Xtrain)
    O2_test = Progeny_pw_act(Xtest)
    
    omics_train = cbind(O1,O2)
    omics_test = cbind(O1_test,O2_test)
    
    colnames(omics_train) = 1:ncol(omics_train)
    colnames(omics_test) = 1:ncol(omics_test)
    
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)))
    Omics = list(omics_train,index,omics_test)
    
  }else if (prod(selected_features == c("Landmark_genes","TF_DoRothEA","Tissue_types"))){
    l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
    O1 = Xtrain[,colnames(Xtrain)%in%l1000_genes]
    O2 = readRDS("Processed_data/S14/DoRothEA_TF.rds")
    O3 = readRDS("Processed_data/S19/sample_tissue_types.rds")
    
    omics_train = cbind(O1,O2,O3)
    colnames(omics_train) = 1:ncol(omics_train)
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)),rep(3,ncol(O3)))
    Omics = list(omics_train,index)
    
  }else if (prod(selected_features == c("Landmark_genes","TF_DoRothEA","TF_decoupleR"))){
    
    l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
    I_G = intersect(l1000_genes,colnames(Xtest))
    
    O1 = Xtrain[,I_G]
    O1_test = Xtest[,I_G]
    
    O2 = DoRothEA(Xtrain)
    O2_test = DoRothEA(Xtest) 
    
    O3 = decoupleR(X = Xtrain, method = "gsva")
    O3_test = decoupleR(X = Xtest, method = "gsva")
    
    omics_train = cbind(O1,O2,O3)
    omics_test = cbind(O1_test,O2_test,O3_test)
    
    colnames(omics_train) = 1:ncol(omics_train)
    colnames(omics_test) = 1:ncol(omics_test)
    
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)),rep(3,ncol(O3)))
    Omics = list(omics_train,index,omics_test)
    
  }else if (prod(selected_features == c("Landmark_genes","Tissue_types","TF_decoupleR"))){
    l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
    O1 = Xtrain[,colnames(Xtrain)%in%l1000_genes]
    O2 = readRDS("Processed_data/S19/sample_tissue_types.rds")
    O3 = decoupleR(X = Xtrain, method = "gsva")
    
    omics_train = cbind(O1,O2,O3)
    colnames(omics_train) = 1:ncol(omics_train)
    
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)),rep(3,ncol(O3)))
    Omics = list(omics_train,index)    
    
  }else if (prod(selected_features == c("Landmark_genes","Tissue_types","TF_DoRothEA","TF_decoupleR"))){
    
    l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
    O1 = Xtrain[,colnames(Xtrain)%in%l1000_genes]
    O2 = readRDS("Processed_data/S19/sample_tissue_types.rds")
    O3 = readRDS("Processed_data/S14/DoRothEA_TF.rds")
    O4 = decoupleR(X = Xtrain, method = "gsva")
    
    omics_train = cbind(O1,O2,O3,O4)
    colnames(omics_train) = 1:ncol(omics_train)
    index = c(rep(1,ncol(O1)),rep(2,ncol(O2)),rep(3,ncol(O3)),rep(4,ncol(O4)))
    Omics = list(omics_train,index)    
  }
  return(Omics)
}
