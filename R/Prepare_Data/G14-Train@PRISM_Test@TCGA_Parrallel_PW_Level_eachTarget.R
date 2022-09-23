rm(list=ls())

library(parallel)
no_cores = detectCores()
cl = makeCluster(no_cores-2)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")

low_sample_drugs = c(1,3,4,5,6,7,12,13,14,19,20,21,22,23,24,25,28,29,31,33,36,37
                     ,39,40,41,42,43,45,47,49,51,52,56,57,58)
sen_PRISM = readRDS("Processed_data/S23/sensitivity_matrix_PRISM_with@TCGA@drugs.rds")
res_TCGA = readRDS("Processed_data/S24/Drug_response_TCGA_binarized.rds")
#res_TCGA = readRDS("Processed_data/S23/Drug_response_matrix_TCGA.rds")
sen_PRISM = sen_PRISM[,-low_sample_drugs]
res_TCGA = res_TCGA[,-low_sample_drugs]

GE = readRDS("Processed_data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")
GE_TCGA = readRDS("Processed_data/S23/expresion_matrix_TCGA.rds")

# Remove genes whose Q3 is zero
q3_genes = apply(GE_TCGA,2,quantile,prob=0.75)
GE_TCGA = GE_TCGA[,-which(q3_genes==0)]
GE = GE[,-which(q3_genes==0)]

N_drug = ncol(sen_PRISM)
drugs = data.frame(colnames(sen_PRISM))
#saveRDS(drugs,"Processed_data/Other/24_drugs.rds")
d = 3
drug = drugs[d,1]

clusterExport(cl, c("GE","GE_TCGA","sen_PRISM","res_TCGA",
                    "drug","d"))
clusterEvalQ(cl, c(source("F7-RandomForest.R"),
                   source("F6-ENet.R"),source("F8-MLP.R"),
                   source("F10-Ridge.R"),source("F11-SGL.R"),
                   source("F13-Lasso.R"),
                   source("F16-Zscore_Normalization.R"),
                   source("F17-Rank_Normalization.R"),
                   source("F18-Combat_Normalization.R")))

LevelLoop = function(i){
  
  print(paste0("The level number is: ", as.character(i)))
  
  #Drug Pathway feature selection
  source("F22-Drug_Pathway_Level_genes_eachTarget.R")
  pathway_gene_set = Drug_Pathway_gene_set_eachTarget(drug = drug, level=i)
  
  if(!isEmpty(pathway_gene_set)){
    
    Result = c()
    for(j in names(pathway_gene_set)){
      
      if(!is.null(pathway_gene_set[[j]])){
        I = intersect(colnames(GE),pathway_gene_set[[j]])
        X = GE[,I]
        X_TCGA = GE_TCGA[,I]
        if(!is.null(ncol(X))){
          index = rep(1,ncol(X))
          ##################
          #X = GE
          #X_TCGA = GE_TCGA
          Xtrain = X[!is.na(sen_PRISM[,d]),]
          ytrain = sen_PRISM[!is.na(sen_PRISM[,d]),d]
          
          Xtest = X_TCGA[!is.na(res_TCGA[,d]),]
          ytest = res_TCGA[!is.na(res_TCGA[,d]),d]
          
          length(ytest)
          if(length(ytest)>10){
            
            #X_Normalization = Rank(Xtrain,Xtest)
            #X_Normalization = Rank(Xtrain,Xtest)
            X_Normalization = Combat_Scale(Xtrain,Xtest)
            
            Xtrain = X_Normalization[[1]]
            Xtest = X_Normalization[[2]]
            N_genes = ncol(Xtrain)
            #source("F15-Feature_Selection_PRISM@TCGA.R")
            #selected_features = c("TF_decoupleR","progeny")
            #Omics_List = Feature_Selection(selected_features,GE = Xtrain ,GE_test = Xtest)
            #Xtrain = Omics_List[[1]]
            #index = Omics_List[[2]]
            #Xtest = Omics_List[[3]]
            
            # Ytrain normalization
            # Mean_ytrain = mean(ytrain)
            # STD_ytrain = sd(ytrain)
            # ytrain = (ytrain-Mean_ytrain)/STD_ytrain
            
            # Models
            #y_pred_Ridge = My_SGL(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest,index = index)
            #y_pred_Ridge = RandomForest(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
            #y_pred_Ridge = ElasticNet(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
            #y_pred_Ridge = Lasso(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
            y_pred_Ridge = Ridge(ytrain = ytrain ,Xtrain = Xtrain, Xtest = Xtest)
            #y_pred_Ridge = MLP(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
            
            # Evaluation
            corr_Ridge = cor(ytest,y_pred_Ridge)
            #print(corr_Ridge)
            #corr_RF = cor(ytest,y_pred_RF)
            #corr_ENet = cor(ytest,y_pred_ENet)
            #corr_Lasso = cor(ytest,y_pred_Lasso)
            #corr_Ridge = cor(ytest , y_pred_Ridge)
            #corr_Ridge = cor(ytest , y_pred_Ridge)
            
            ttest = t.test(y_pred_Ridge[ytest==1], y_pred_Ridge[ytest==2], alternative="greater")$p.value
            Ranksum = wilcox.test(y_pred_Ridge[ytest==1], y_pred_Ridge[ytest==2], alternative ="greater")$p.value
            
        }else{
          N_genes = 0
          corr_Ridge = 0
          ttest = 1
          Ranksum = 1
        }
          
          Result = rbind(Result, cbind(corr_Ridge,ttest,Ranksum,N_genes))
          }else{
            N_genes = 0
            corr_Ridge = 0
            ttest = 1
            Ranksum = 1
            Result = rbind(Result, cbind(corr_Ridge,ttest,Ranksum,N_genes))
          }
      }else{
        Result = c()
      }
    }
      }else{
    #corr_Ridge = 0
    #ttest = 1
    #Ranksum = 1
    Result = c()
  }
  result = cbind(Result, Targets = names(pathway_gene_set))

  return(result)
}
N_Level = 10
result = parLapply(cl, sapply(1:N_Level, list), LevelLoop)

Result = list()
for (k in 1:N_Level){
  if(ncol(result[[k]])!=1)
    Result[[k]] = result[[k]]
}

stopCluster(cl)

df = c()
for(i in 1:length(Result)){
  target_i = data.frame(Result[[i]])
  df = rbind(df,cbind(p_val = round(-log10(as.numeric(target_i$Ranksum)),2), 
                      level = i, target = target_i$Targets,
                      N_genes = target_i$N_genes))
}

df = data.frame(df)
df$p_val = as.numeric(df$p_val)
df$level = as.numeric(df$level)

# plt = ggplot(df) +
#   geom_line(aes(x = level, y = p_val, color = target), size = .5) + 
#   geom_point(aes(x = level, y = p_val, color = target), size = 1) +
#   theme_classic(base_size = 10) + theme(legend.position = "right") +
#   geom_hline(yintercept = -log10(0.05),linetype = "dashed")+
#   labs(title= drug, x ="Level", y = "-log(P-value)")+
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 
# pdf(paste0("Figures/FS/Drug_pathway/", drug,".pdf"), height = 4, width = 5)
# plt
# dev.off()


