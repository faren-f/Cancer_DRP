rm(list=ls())

library(parallel)
no_cores = detectCores()
cl = makeCluster(no_cores-2)

setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")
sen = readRDS("Processed_data/Other/Sen_PRISM_24_Drugs.rds")
GE = readRDS("Processed_data/S23/expresion_matrix_PRISM_with@TCGA@genes.rds")


drugs = data.frame(colnames(sen))
d = 12
drug = drugs[d,1]

clusterExport(cl, c("GE","sen","drug","d"))
clusterEvalQ(cl, c(library(caTools), source("F7-RandomForest.R"),
                   source("F6-ENet.R"),source("F8-MLP.R"),
                   source("F10-Ridge.R"),source("F11-SGL.R"),
                   source("F13-Lasso.R"),
                   source("F16-Zscore_Normalization.R"),
                   source("F17-Rank_Normalization.R")))

LevelLoop = function(i){
  
  print(paste0("The level number is: ", as.character(i)))
  
  source("F22-Drug_Pathway_Level_genes_eachTarget.R")
  pathway_gene_set = Drug_Pathway_gene_set_eachTarget(drug = drug, level=i)
  
  if(!isEmpty(pathway_gene_set)){
    
    Result = c()
    for(j in names(pathway_gene_set)){
      
      if(!is.null(pathway_gene_set[[j]])){
        I = intersect(colnames(GE),pathway_gene_set[[j]])
        X = GE[,I]
        
        if(!is.null(ncol(X))){
          index = rep(1,ncol(X))
          ##################
          
          X = X[!is.na(sen[,d]),]
          y = sen[!is.na(sen[,d]),d]
 
          
          if(length(ytest)>10){
            
            sample = sample.split(y, SplitRatio = .8)
            
            Xtrain = subset(X, sample == TRUE)
            Xtest  = subset(X, sample == FALSE)
            ytrain = subset(y, sample == TRUE)
            ytest  = subset(y, sample == FALSE)
            
    
            # Models
            #y_pred_Ridge = My_SGL(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest,index = index)
            #y_pred_Ridge = RandomForest(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
            #y_pred_Ridge = ElasticNet(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
            #y_pred_Ridge = Lasso(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
            #y_pred_Ridge = Ridge(ytrain = ytrain ,Xtrain = Xtrain, Xtest = Xtest)
            #y_pred_Ridge = MLP(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
            
            y_pred_Ridge = Ridge(ytrain = ytrain ,Xtrain = Xtrain, 
                                   Xtest = Xtest)
            
            # y_pred_Ridge = MLP(ytrain = ytrain ,Xtrain = Xtrain, 
            #                        Xtest = Xtest)
            
            # Evaluation
            corr_Ridge = cor(ytest,y_pred_Ridge)
            corr_Ridge = cor(ytest,y_pred_Ridge)
            plot(ytest,y_pred_Ridge,ylim = c(0,1.8),xlim = c(0,1.8))
            plot(ytest,y_pred_Ridge)
            
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


