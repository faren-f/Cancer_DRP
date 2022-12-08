rm(list = ls())
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
sen = readRDS("Processed_data/S1/sensitivity_matrix_AUC.rds")
response = read.csv("Raw_data/PRISM/Secondary/secondary-screen-dose-response-curve-parameters.csv")

N_drugs = 1448
RF_Landmark = c()
ENet_Landmark = c()
Lasso_Landmark = c()
Ridge_Landmark = c()
MLP_Landmark = c()

for(i in 1:N_drugs){
  R_Landmark = readRDS(paste0("Processed_from_SLURM/Results_Landmark_All_Models/Result_",as.character(i),".rds"))
  
  Ridge_Landmark = rbind(Ridge_Landmark, R_Landmark[4,])
  MLP_Landmark = rbind(MLP_Landmark, R_Landmark[5,])
  Lasso_Landmark = rbind(Lasso_Landmark, R_Landmark[3,])
  ENet_Landmark = rbind(ENet_Landmark, R_Landmark[2,])
  RF_Landmark = rbind(RF_Landmark, R_Landmark[1,])
}

rownames(Ridge_Landmark) = colnames(sen)
rownames(MLP_Landmark) = colnames(sen)
rownames(Lasso_Landmark) = colnames(sen)
rownames(ENet_Landmark) = colnames(sen)
rownames(RF_Landmark) = colnames(sen)

# Finding drugs with the same mechanism of actions
Moa = unique(response$moa)
Moa = strsplit(Moa,",")
Moa = unique(unlist(Moa))

w = strsplit(response$moa, ",")
I_NA = which(is.na(w))
for(r in I_NA){
  w[[r]] = "Unknown"
}

d = list()
l = c()
D = c()
for(m in Moa){
  d[m] = list(unique(response[which(sapply(w, FUN=function(x) m %in% x)),12]))
  l = c(l, length(d[[m]]))
  if(length(d[[m]])>5){
    D = c(D, m)
  }
}
#saveRDS(d, "Final_Result/List_of_Drug_MOA.rds")

boxplot(Ridge_Landmark[d[[D[1]]],1], Ridge_Landmark[d[[D[2]]],1],
        Ridge_Landmark[d[[D[3]]],1], Ridge_Landmark[d[[D[4]]],1],
        Ridge_Landmark[d[[D[5]]],1], Ridge_Landmark[d[[D[6]]],1],
        Ridge_Landmark[d[[D[7]]],1], Ridge_Landmark[d[[D[8]]],1],
        Ridge_Landmark[d[[D[9]]],1], Ridge_Landmark[d[[D[10]]],1],
        Ridge_Landmark[d[[D[11]]],1], Ridge_Landmark[d[[D[12]]],1],
        Ridge_Landmark[d[[D[13]]],1], Ridge_Landmark[d[[D[14]]],1],
        Ridge_Landmark[d[[D[15]]],1], Ridge_Landmark[d[[D[16]]],1],
        Ridge_Landmark[d[[D[17]]],1], Ridge_Landmark[d[[D[18]]],1],
        Ridge_Landmark[d[[D[19]]],1], Ridge_Landmark[d[[D[20]]],1],
        names=NA, cex=.5)
        

boxplot(Ridge_Landmark[d[[D[21]]],1], Ridge_Landmark[d[[D[22]]],1],
        Ridge_Landmark[d[[D[23]]],1], Ridge_Landmark[d[[D[24]]],1],
        Ridge_Landmark[d[[D[25]]],1], Ridge_Landmark[d[[D[26]]],1],
        Ridge_Landmark[d[[D[27]]],1], Ridge_Landmark[d[[D[28]]],1],
        Ridge_Landmark[d[[D[29]]],1], Ridge_Landmark[d[[D[30]]],1],
        Ridge_Landmark[d[[D[31]]],1], Ridge_Landmark[d[[D[32]]],1],
        Ridge_Landmark[d[[D[33]]],1], Ridge_Landmark[d[[D[34]]],1],
        Ridge_Landmark[d[[D[35]]],1], Ridge_Landmark[d[[D[36]]],1],
        Ridge_Landmark[d[[D[37]]],1], Ridge_Landmark[d[[D[38]]],1],
        Ridge_Landmark[d[[D[39]]],1], Ridge_Landmark[d[[D[40]]],1],
        names=NA, cex=.5)

boxplot(Ridge_Landmark[d[[D[41]]],1], Ridge_Landmark[d[[D[42]]],1],
        Ridge_Landmark[d[[D[43]]],1], Ridge_Landmark[d[[D[44]]],1],
        Ridge_Landmark[d[[D[45]]],1], Ridge_Landmark[d[[D[46]]],1],
        Ridge_Landmark[d[[D[47]]],1], Ridge_Landmark[d[[D[48]]],1],
        Ridge_Landmark[d[[D[49]]],1], Ridge_Landmark[d[[D[40]]],1],
        Ridge_Landmark[d[[D[51]]],1], Ridge_Landmark[d[[D[52]]],1],
        Ridge_Landmark[d[[D[53]]],1], Ridge_Landmark[d[[D[54]]],1],
        Ridge_Landmark[d[[D[55]]],1], Ridge_Landmark[d[[D[56]]],1],
        Ridge_Landmark[d[[D[57]]],1], Ridge_Landmark[d[[D[58]]],1],
        Ridge_Landmark[d[[D[59]]],1], Ridge_Landmark[d[[D[60]]],1],
        names=NA, cex=.5)

boxplot(Ridge_Landmark[d[[D[61]]],1], Ridge_Landmark[d[[D[62]]],1],
        Ridge_Landmark[d[[D[63]]],1], Ridge_Landmark[d[[D[64]]],1],
        Ridge_Landmark[d[[D[65]]],1], Ridge_Landmark[d[[D[66]]],1],
        Ridge_Landmark[d[[D[67]]],1], Ridge_Landmark[d[[D[68]]],1],
        Ridge_Landmark[d[[D[69]]],1], Ridge_Landmark[d[[D[70]]],1],
        Ridge_Landmark[d[[D[71]]],1], Ridge_Landmark[d[[D[72]]],1],
        Ridge_Landmark[d[[D[73]]],1], Ridge_Landmark[d[[D[74]]],1],
        Ridge_Landmark[d[[D[75]]],1], Ridge_Landmark[d[[D[76]]],1],
        Ridge_Landmark[d[[D[77]]],1], Ridge_Landmark[d[[D[78]]],1],
        Ridge_Landmark[d[[D[79]]],1], Ridge_Landmark[d[[D[80]]],1],
        names=NA, cex=.5)


boxplot(Ridge_Landmark[d[[D[81]]],1], Ridge_Landmark[d[[D[82]]],1],
        Ridge_Landmark[d[[D[83]]],1], Ridge_Landmark[d[[D[84]]],1],
        Ridge_Landmark[d[[D[85]]],1], Ridge_Landmark[d[[D[86]]],1],
        Ridge_Landmark[d[[D[87]]],1], Ridge_Landmark[d[[D[88]]],1],
        Ridge_Landmark[d[[D[89]]],1], Ridge_Landmark[d[[D[90]]],1],
        Ridge_Landmark[d[[D[91]]],1], Ridge_Landmark[d[[D[92]]],1],
        Ridge_Landmark[d[[D[93]]],1], Ridge_Landmark[d[[D[94]]],1],
        Ridge_Landmark[d[[D[95]]],1], Ridge_Landmark[d[[D[96]]],1],
        Ridge_Landmark[d[[D[97]]],1], Ridge_Landmark[d[[D[98]]],1],
        Ridge_Landmark[d[[D[99]]],1], Ridge_Landmark[d[[D[100]]],1],
        names=NA, cex=.5)

boxplot(Ridge_Landmark[d[[D[101]]],1], Ridge_Landmark[d[[D[102]]],1],
        Ridge_Landmark[d[[D[103]]],1], Ridge_Landmark[d[[D[104]]],1],
        Ridge_Landmark[d[[D[105]]],1], Ridge_Landmark[d[[D[106]]],1],
        Ridge_Landmark[d[[D[107]]],1], Ridge_Landmark[d[[D[108]]],1],
        Ridge_Landmark[d[[D[109]]],1], Ridge_Landmark[d[[D[110]]],1],
        Ridge_Landmark[d[[D[111]]],1], Ridge_Landmark[d[[D[112]]],1],
        Ridge_Landmark[d[[D[113]]],1], Ridge_Landmark[d[[D[114]]],1],
        Ridge_Landmark[d[[D[115]]],1], Ridge_Landmark[d[[D[116]]],1],
        Ridge_Landmark[d[[D[117]]],1], Ridge_Landmark[d[[D[118]]],1],
        Ridge_Landmark[d[[D[119]]],1], Ridge_Landmark[d[[D[120]]],1],
        names=NA, cex=.5)

boxplot(Ridge_Landmark[d[[D[121]]],1], Ridge_Landmark[d[[D[122]]],1],
        Ridge_Landmark[d[[D[123]]],1], Ridge_Landmark[d[[D[124]]],1],
        Ridge_Landmark[d[[D[125]]],1], Ridge_Landmark[d[[D[126]]],1],
        Ridge_Landmark[d[[D[127]]],1], Ridge_Landmark[d[[D[128]]],1],
        Ridge_Landmark[d[[D[129]]],1], Ridge_Landmark[d[[D[130]]],1],
        Ridge_Landmark[d[[D[131]]],1], Ridge_Landmark[d[[D[132]]],1],
        Ridge_Landmark[d[[D[133]]],1], Ridge_Landmark[d[[D[134]]],1],
        Ridge_Landmark[d[[D[135]]],1], Ridge_Landmark[d[[D[136]]],1],
        Ridge_Landmark[d[[D[137]]],1], Ridge_Landmark[d[[D[138]]],1],
        Ridge_Landmark[d[[D[139]]],1], Ridge_Landmark[d[[D[140]]],1],
        names=NA, cex=.5)





boxplot(Ridge_Landmark[d[[141]],1], Ridge_Landmark[d[[142]],1],
        Ridge_Landmark[d[[143]],1], Ridge_Landmark[d[[144]],1],
        Ridge_Landmark[d[[145]],1], Ridge_Landmark[d[[146]],1],
        Ridge_Landmark[d[[147]],1], Ridge_Landmark[d[[148]],1],
        Ridge_Landmark[d[[149]],1], Ridge_Landmark[d[[150]],1],
        Ridge_Landmark[d[[151]],1], Ridge_Landmark[d[[152]],1],
        Ridge_Landmark[d[[153]],1], Ridge_Landmark[d[[154]],1],
        Ridge_Landmark[d[[155]],1], Ridge_Landmark[d[[156]],1],
        Ridge_Landmark[d[[157]],1], Ridge_Landmark[d[[158]],1],
        Ridge_Landmark[d[[159]],1], Ridge_Landmark[d[[160]],1])

boxplot(Ridge_Landmark[d[[161]],1], Ridge_Landmark[d[[162]],1],
        Ridge_Landmark[d[[163]],1], Ridge_Landmark[d[[164]],1],
        Ridge_Landmark[d[[165]],1], Ridge_Landmark[d[[166]],1],
        Ridge_Landmark[d[[167]],1], Ridge_Landmark[d[[168]],1],
        Ridge_Landmark[d[[169]],1], Ridge_Landmark[d[[170]],1],
        Ridge_Landmark[d[[171]],1], Ridge_Landmark[d[[172]],1],
        Ridge_Landmark[d[[173]],1], Ridge_Landmark[d[[174]],1],
        Ridge_Landmark[d[[175]],1], Ridge_Landmark[d[[176]],1],
        Ridge_Landmark[d[[177]],1], Ridge_Landmark[d[[178]],1],
        Ridge_Landmark[d[[179]],1], Ridge_Landmark[d[[180]],1])

boxplot(Ridge_Landmark[d[[1]],1], Ridge_Landmark[d[[2]],1],
        Ridge_Landmark[d[[3]],1], Ridge_Landmark[d[[4]],1],
        Ridge_Landmark[d[[5]],1], Ridge_Landmark[d[[6]],1],
        Ridge_Landmark[d[[7]],1], Ridge_Landmark[d[[8]],1],
        Ridge_Landmark[d[[9]],1], Ridge_Landmark[d[[0]],1],
        Ridge_Landmark[d[[1]],1], Ridge_Landmark[d[[2]],1],
        Ridge_Landmark[d[[3]],1], Ridge_Landmark[d[[4]],1],
        Ridge_Landmark[d[[5]],1], Ridge_Landmark[d[[6]],1],
        Ridge_Landmark[d[[7]],1], Ridge_Landmark[d[[8]],1],
        Ridge_Landmark[d[[9]],1], Ridge_Landmark[d[[0]],1])

boxplot(Ridge_Landmark[d[[1]],1], Ridge_Landmark[d[[2]],1],
        Ridge_Landmark[d[[3]],1], Ridge_Landmark[d[[4]],1],
        Ridge_Landmark[d[[5]],1], Ridge_Landmark[d[[6]],1],
        Ridge_Landmark[d[[7]],1], Ridge_Landmark[d[[8]],1],
        Ridge_Landmark[d[[9]],1], Ridge_Landmark[d[[0]],1],
        Ridge_Landmark[d[[1]],1], Ridge_Landmark[d[[2]],1],
        Ridge_Landmark[d[[3]],1], Ridge_Landmark[d[[4]],1],
        Ridge_Landmark[d[[5]],1], Ridge_Landmark[d[[6]],1],
        Ridge_Landmark[d[[7]],1], Ridge_Landmark[d[[8]],1],
        Ridge_Landmark[d[[9]],1], Ridge_Landmark[d[[0]],1])

boxplot(Ridge_Landmark[d[[1]],1], Ridge_Landmark[d[[2]],1],
        Ridge_Landmark[d[[3]],1], Ridge_Landmark[d[[4]],1],
        Ridge_Landmark[d[[5]],1], Ridge_Landmark[d[[6]],1],
        Ridge_Landmark[d[[7]],1], Ridge_Landmark[d[[8]],1],
        Ridge_Landmark[d[[9]],1], Ridge_Landmark[d[[0]],1],
        Ridge_Landmark[d[[1]],1], Ridge_Landmark[d[[2]],1],
        Ridge_Landmark[d[[3]],1], Ridge_Landmark[d[[4]],1],
        Ridge_Landmark[d[[5]],1], Ridge_Landmark[d[[6]],1],
        Ridge_Landmark[d[[7]],1], Ridge_Landmark[d[[8]],1],
        Ridge_Landmark[d[[9]],1], Ridge_Landmark[d[[0]],1])










#boxplot(Results_Ridge, names=NA, cex=.5)
