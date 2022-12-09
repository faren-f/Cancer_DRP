rm(list=ls())

library(caTools)
source("Functions/F10-Ridge.R")
source("Functions/F6-ENet.R")
source("Functions/F8-MLP.R")
source("Functions/F13-Lasso.R")
source("Functions/F7-RandomForest.R")

### SLURM_ARRAY_TASK_ID
task_id = Sys.getenv("SLURM_ARRAY_TASK_ID")
task_id = as.numeric(task_id)
i = task_id

sen = readRDS("Data/sensitivity_matrix_AUC.rds")
GE = readRDS("Data/expresion_matrix.rds")

#boxplot(GE[,1:20], names=NA, cex=.1, outline = FALSE, main="background corrected data")

l1000_genes = readRDS("Processed_Data/S18/Landmark_genes.rds")
GE = GE[,colnames(GE)%in%l1000_genes]

Models = c("RandomForest","ElasticNet", "Lasso","Ridge","MLP")

#for (i in 1:ncol(sen)){            # drug loop
print(paste0("The drug number is: ", as.character(i)))
Mean_Corr = c()
STD_Corr = c()

Mean_MSE = c()
STD_MSE = c()

for (M in Models){             # model loop
  model = get(M)
  Corr = c()
  MSE = c()
  print(M)
  
  for (j in 1:50){           # repeat loop
    print(paste0("The repeat number is: ", as.character(j)))
    
    X = GE[!is.na(sen[,i]),]
    y = sen[!is.na(sen[,i]),i]
    
    # normalization
    X = scale(X)
    y = scale(y)
    y = y[,1]
    
    sample = sample.split(y, SplitRatio = .8)
    
    Xtrain = subset(X, sample == TRUE)
    Xtest  = subset(X, sample == FALSE)
    ytrain = subset(y, sample == TRUE)
    ytest  = subset(y, sample == FALSE)
    
    # Models
    y_pred = model(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    
    # Evaluation
    corr = cor(ytest,y_pred)
    Corr = c(Corr, corr)
    
    mse =  mean((ytest - y_pred)^2)
    MSE = c(MSE, mse)
  }
  
  Mean_Corr = c(Mean_Corr, mean(Corr))
  STD_Corr = c(STD_Corr, sd(Corr))
  
  Mean_MSE = c(Mean_MSE, mean(MSE))
  STD_MSE = c(STD_MSE, sd(MSE))
  
}

Result = data.frame(Mean_Corr = Mean_Corr, STD_Corr = STD_Corr,
                    Mean_MSE = Mean_MSE, STD_MSE = STD_MSE)
rownames(Result) = Models
  
saveRDS(Result, paste0("Results/temp/Result_",as.character(i),".rds"))


