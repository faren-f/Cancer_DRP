rm(list=ls())

Results = readRDS("All_Results/Result_RF_All_Drugs.rds")
good_Results = data.frame(index = which(Results[,1]>0.3))
saveRDS(good_Results,"All_Results/good_drugs_in_PRISM.rds")

sen = readRDS("Processed_Data/S1/sensitivity_matrix.rds")
sen_re = sen[,good_Results[,1]]
saveRDS(sen_re,"All_Results/sen_PRISM_good_drugs.rds")
