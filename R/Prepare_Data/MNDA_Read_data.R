rm(list=ls())

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
response = read.csv("Raw_data/PRISM/Secondary/secondary-screen-dose-response-curve-parameters.csv")
sen_PRISM = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")
GE_PRISM = readRDS("Processed_data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")
q3_genes = apply(GE_PRISM,2,quantile,prob=0.75)
GE_PRISM = GE_PRISM[,-which(q3_genes==0)]

tissues = strsplit(response$ccle_name, "_")
response[,21] = sapply(tissues,FUN = function(x){return(x[2])})

#drug = response[response$name == "anastrozole",]
#drug_tissue = drug[drug[,21] == "LUNG",]
drug = response[response$name == "tamoxifen",]
drug_tissue = drug[drug[,21] == "LUNG",]
drug_tissue = drug_tissue[!is.na(drug_tissue$depmap_id),]
cell_id = unique(drug_tissue[,2])

AUC = matrix(0,length(cell_id),1)
rownames(AUC) = cell_id
c = 0
for (i in cell_id){
  c = c+1
  print(c)
  
  cell_i = which(drug_tissue$depmap_id == i)
    
  if (length(cell_i) == 1){
    AUC[i,1] = drug_tissue[cell_i,"auc"]
    
  }else{
    AUC[i,1] = mean(drug_tissue[cell_i,"auc"])
    
  }
}

hist(AUC[,1],30)

plot(sort(AUC),type = "l")
abline(v = 36, col = "red")
abline(h = 1.15, col = "red")
abline(h = 1.25, col = "red")

sorted_auc = sort(AUC)
a = 2
dif = c()
for(i in seq(1,(length(sorted_auc)-a),a)){
  dif = c(dif, sorted_auc[i+a]-sorted_auc[i])
}

plot(dif,type = "l")
a*which.min(dif)

intesected_celllines = intersect(rownames(GE_PRISM),rownames(AUC))
AUC = data.frame(AUC[intesected_celllines,])

res = rownames(AUC)[which(AUC<1.15)]
non_res = rownames(AUC)[which(AUC>1.25)]


y = c(rep("res",length(res)),rep("non_res",length(non_res)))
GE_res = GE_PRISM[rownames(GE_PRISM) %in% res,]
GE_non_res = GE_PRISM[rownames(GE_PRISM) %in% non_res,]
GE = rbind(GE_res,GE_non_res)
Var = apply(GE,2,var)
hist(Var,100)
a = order(Var, decreasing = TRUE)
GE = GE[,a[1:2000]]
Var = apply(GE,2,var)
hist(Var,100)

GE = scale(GE)

data = list(Expr = GE, response = y)
saveRDS(data,"MNDA/data_lung_tamoxifen_2000genes.rds")

