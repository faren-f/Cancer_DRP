rm(list=ls())
setwd("~/Desktop/Programming/R_Root/PhD_Project/DRP_PRISM/")

# Read_Data ---------------------------------------------------------------
emb_celline = read.csv("Data/Embedding/RGCN/Embedding_cellline_prod_1L.csv", header = F)
emb_drug = read.csv("Data/Embedding/RGCN/Embedding_drug_prod_1L.csv",header = F)
GE = readRDS("Data/Processed_Data/expresion_matrix.rds")
sen = readRDS("Data/Processed_Data/sensitivity_matrix.rds")
Fingerprints = readRDS("Data/Processed_Data/Fingerprints.rds")  


# cell line similarity based on sensitivity and embedding
I = c()
for (j in 1:nrow(sen)){
  if((sum(!is.na(sen[j,]))>1000) & (sum(!is.na(sen[j,]))<1448)) {
    I = c(I, j)
  }}

#sen_1 = sen[I,] 
#emb_celline_1 = emb_celline[I,]
cellline_similarity_emb = as.matrix(dist(emb_celline))
cellline_similarity_sen = as.matrix(dist(sen))


#### different distance measures
#cellline_similarity_emb = cor(t(emb_celline))
#cellline_similarity_emb = as.matrix(dist(emb_celline))
#cellline_similarity_emb = as.matrix(emb_celline) %*% t(emb_celline)
#cellline_similarity_emb = cellline_similarity_emb / apply(cellline_similarity_emb, 2, function(x){return(sqrt(sum(x^2)))})

#cellline_similarity_sen = cor(t(sen), use="complete.obs")
#cellline_similarity_sen = as.matrix(dist(sen))
#cellline_similarity_sen = as.matrix(sen) %*% t(sen)
####

cor(matrix(cellline_similarity_emb,length(cellline_similarity_emb)),
    matrix(cellline_similarity_sen,length(cellline_similarity_emb)))

plot(cellline_similarity_emb,cellline_similarity_sen, pch = 20, cex = .01)

# cell line similarity based on gene expression and embedding

var_GE = apply(GE,2,var)
var_GE = sort(var_GE,decreasing = TRUE)
hist(var_GE,100)
thr_var_GE = 0.8
abline(v = thr_var_GE, col = "red")

GE = GE[,var_GE > thr_var_GE]

# Normalization Gene Expresion
GE = scale(GE)
cellline_similarity_emb = as.matrix(dist(emb_celline))
cellline_similarity_GE = as.matrix(dist(GE))

cor(matrix(cellline_similarity_emb,length(cellline_similarity_emb)),
matrix(cellline_similarity_GE,length(cellline_similarity_emb)))

plot(cellline_similarity_emb,cellline_similarity_GE, pch = 20, cex = .1)

## Drug similarity based on fingerprints and embedding
FP = matrix(0, length(Fingerprints), 1024)
for (i in 1:length(Fingerprints)){
  FP_bits_on = Fingerprints[[i]]
  FP[i,FP_bits_on@bits] = 1
}

drug_similarity_emb = as.matrix(dist(emb_drug))
drug_similarity_FP = as.matrix(dist(FP))

cor(matrix(drug_similarity_emb,length(drug_similarity_emb)),
    matrix(drug_similarity_FP,length(drug_similarity_emb)))

plot(cellline_similarity_emb,drug_similarity_FP, pch = 20, cex = .1)

# drug similarity based on sensitivity and embedding
drug_similarity_emb = as.matrix(dist(emb_drug))
drug_similarity_sen = as.matrix(dist(t(sen)))

cor(matrix(drug_similarity_emb,length(drug_similarity_emb)),
    matrix(drug_similarity_sen,length(drug_similarity_emb)))

plot(drug_similarity_emb,drug_similarity_sen, pch = 20, cex = .1)

# drug similarity based on sensitivity and Fingerprint
drug_similarity_sen = as.matrix(dist(t(sen)))

cor(matrix(drug_similarity_FP,length(drug_similarity_FP)),
    matrix(drug_similarity_sen,length(drug_similarity_FP)))

plot(drug_similarity_FP,drug_similarity_sen, pch = 20, cex = .1)



# Different Similarity Calculation
cellline_similarity_emb = cor(t(emb_celline))
cellline_similarity_emb = as.matrix(dist(emb_celline))
cellline_similarity_emb = as.matrix(emb_celline) %*% t(emb_celline)
cellline_similarity_emb = cellline_similarity_emb / apply(cellline_similarity_emb, 2, function(x){return(sqrt(sum(x^2)))})

cellline_similarity_sen = cor(t(sen), use="complete.obs")
cellline_similarity_sen = as.matrix(dist(sen))
cellline_similarity_sen = as.matrix(sen) %*% t(sen)

#cellline_similarity_sen[cellline_similarity_emb>10] = 0
#cellline_similarity_emb[cellline_similarity_emb>10] = 0

cor(matrix(cellline_similarity_emb,length(cellline_similarity_emb)),
    matrix(cellline_similarity_sen,length(cellline_similarity_emb)))
plot(cellline_similarity_emb,cellline_similarity_sen, pch = 20, cex = .1,xlim = c(2,12),ylim = c(2,12))


y = as.matrix(emb_celline) %*% t(emb_drug)
cor(matrix(y,length(y)), matrix(sen,length(y)), use="complete.obs")
plot(y,sen, pch = 20, cex = .01)















#PCA = prcomp(t(emb_celline), scale. = TRUE)
#Qc = PCA$rotation[,1:2]
#plot(Qc, cex = .5, pch = 20)

#PCA = prcomp(t(emb_drug), scale. = TRUE)
#Qd = PCA$rotation[,1:2]
#plot(Qd, cex = .5, pch = 20)

#plot(Qc, cex = .3, pch = 20, col = "red", xlim = c(.02,.06), ylim = c(-.06,.1))
#par(new =  TRUE)
#plot(Qd, cex = .2, pch = 20, col = "blue", xlim = c(.02,.06),  ylim = c(-.06,.1))

