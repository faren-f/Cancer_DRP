rm(list=ls())
setwd("~/Desktop/Codes/Cancer_DRP/R/")


# Library -----------------------------------------------------------------

library(webchem)
library(rcdk)

# Read_Data ---------------------------------------------------------------

expr = readRDS("Data/Processed_Data/expresion_matrix.rds")
expr_norm = readRDS("Data/Processed_Data/expresion_normalized_matrix.rds")

sen = readRDS("Data/Processed_Data/sensitivity_matrix.rds")
response = read.csv("Data/PRISM_Raw_Dataset/Secondary/secondary-screen-dose-response-curve-parameters.csv")


# Similarity calculation --------------------------------------------------


# Cellline-Cellline.Similarity --------------------------------------------
#'@Cellline-Cellline.Similarity

#1)cell line-cell line similarity across genes
sample_sim_exp = cor(t(expr))
#saveRDS(sample_sim_exp,"Processed_Data/Sample_Sim_Exp.rds")

#2)cell line-cell line similarity across drugs
sample_sim_sen = cor(t(sen), use="complete.obs") # "use" is ignoring "NA"

#3)normalized cell line-cell line similarity across genes
#sample_sim_exp_norm = cor(t(expr_norm))

## Plot cell line-cell line similarity (1) against (2)
concat_sample_sim_exp = matrix(sample_sim_exp, length(sample_sim_exp))
concat_sample_sim_sen = matrix(sample_sim_sen, length(sample_sim_sen))

#pdf(file = "Figures/Fig3_A_Sample_Sim.pdf", width = 7, height = 7)
plot(concat_sample_sim_exp, concat_sample_sim_sen, cex=.005, pch=20,
     xlab = "Expresion-based similarity", 
     ylab = "Drug response-based similarity")
#dev.off()

# Drug-Drug.Similarity ----------------------------------------------------

#'@Drug-Drug.Similarity

# CID_Name ----------------------------------------------------------------
### obtain CID numbers of drugs using their names
#drug_name = unique(response$name)
#CID_Name_all = get_cid(drug_name, from = "name",
             # domain = c("compound", "substance", "assay"),
             # match = c("all", "first", "ask", "na"),
             # verbose = getOption("verbose"),arg = NULL, first = NULL)

#CID_Name_all = cbind(CID_Name_all[,2],CID_Name_all[,1])
#saveRDS(CID_Name_all, "Processed_Data/CID_Name_all.rds")
#CID_Name_all = readRDS("Processed_Data/CID_Name_all.rds")
#CID_Name = CID_Name_all[!duplicated(CID_Name_all[,1]),]
#CID_Name = data.frame(CID_Name)
### There are some "NA" in the compound_CIDs, so they are hand curated using 
                  ###Pubchem and ChEMBL datasets.   
# CID_Name[512,2] = "46843057"
# CID_Name[600,2] = "2249"
# CID_Name[674,2] = "260439"
# CID_Name[717,2] = "512282"
# CID_Name[877,2] = "6603754"
# CID_Name[1028,2] = "163838"
# CID_Name[1128,2] = "672296"
# CID_Name[1227,2] = "2210370"
# CID_Name[1381,2] = "461310"
# CID_Name[213,2] = "2724387"
# CID_Name[615,2] = "16051930"
# which(is.na(CID_Name[,2]))
#CID_Name = cbind(CID_Name[,2],CID_Name[,1])
#colnames(CID_Name) = c("CID","Name")
#saveRDS(CID_Name,"Processed_Data/CID_Name.rds")
CID_Name = readRDS("Data/Processed_Data/CID_Name.rds")

# CID_Smile_PubChem -------------------------------------------------------
######CID_Smile is downloaded from https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/
###1)Smiles from PubChem: It did not work
#CID_Smile_all = read.table("Drug_Feature_Data/CID-SMILES", sep = "\t")
#saveRDS(CID_Smile_all, "Data/Processed_Data/CID_Smile_all.rds")

## I saved CID_Smile_all with the "rds" format since it is faster to read in R 
#CID_Smile_all = readRDS("Processed_Data/CID_Smile_all.rds")
#length(intersect(CID_Smile[,1],as.numeric(CID_Name[,2])))              
#CID_Smile = CID_Smile_all[CID_Smile_all[,1] %in% as.numeric(CID_Name[,2]),]
#saveRDS(CID_Smile, "Data/Processed_Data/CID_Smile.rds")
#CID_Smile = readRDS("Data/Processed_Data/CID_Smile.rds")

#Smiles =  CID_Smile[,2]
#molecule <- parse.smiles(Smiles)
## There are 52 warnings. maybe some Smiles are not correct 
## Warned Indeces 
#ind_warning = c()
#for (i in 1:length(molecule))
  #if (length(molecule[[i]])==0)
    #ind_warning = c(ind_warning, i)

#CID_Smile_warned = CID_Smile[ind_warning,]
######## modify warned smiles
#CID_Smile[32,2] = "CC(C)(C#N)C1=CC(=CC(=C1)CN2C=NC=N2)C(C)(C)C#N"
#CID_Smile[42,2] = "CC(O)(CS(=O)(=O)c1ccc(F)cc1)C(=O)Nc1ccc(C#N)c(C(F)(F)F)c1"
#CID_Smile[104,2] = "CCOC(=O)[C@@H](C#N)[C@@H]1c2cc(Br)ccc2OC(N)=C1C(=O)OCC"
## Instead of modifying these warnings we try to use Smiles from PRISM


# CID_Smile_PRISM ---------------------------------------------------------
### 2) Smiles from PRISM: It works
Smile = response$smiles
Smile = sapply(strsplit(Smile, split = ",", perl=T), 
               FUN = function(x){return(x[1])})

Name_Smile = cbind(Name = response$name,Smile = Smile)
Name_Smile = Name_Smile[!duplicated(Name_Smile[,1]),]

### To save CID number, Name and Smiles of drugs in a dataframe 
CID_Name_Smile = merge(Name_Smile,CID_Name,by = "Name")
CID_Name_Smile = CID_Name_Smile[,c(3,1,2)]
rownames(CID_Name_Smile) = CID_Name_Smile[,2]
CID_Name_Smile = CID_Name_Smile[colnames(sen),]


#saveRDS(CID_Name_Smile,"Processed_Data/CID_Name_Smile.rds")

### get fingerprints using "rcdk" package
Smiles =  CID_Name_Smile[,3]
molecule <- parse.smiles(Smiles)
Fingerprints = lapply(molecule, get.fingerprint, type='circular')
saveRDS(Fingerprints,"Data/Processed_Data/Fingerprints.rds")  

Fingerprints_sim = fingerprint::fp.sim.matrix(Fingerprints, method='tanimoto')
rownames (Fingerprints_sim) =  CID_Name_Smile[,2]
colnames (Fingerprints_sim) =  CID_Name_Smile[,2]

saveRDS(Fingerprints_sim,"Data/Processed_Data/Fingerprints_sim.rds")  

## Just for visualization
## convert similarity to distance(disimilarity) for clustring
Fingerprints_dist = 1-Fingerprints_sim
pheatmap::pheatmap(Fingerprints_dist)



