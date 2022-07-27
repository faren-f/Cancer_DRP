rm(list = ls())
library(keras)
library(tensorflow)

library(corrplot)
library(ggvis)
library(caTools)
setwd("~/Desktop/Cancer_DRP/R/Single_Drug/")
## Read data
GE = readRDS("Raw_data/expresion_matrix.rds")
sen = readRDS("Raw_data/sensitivity_matrix.rds")


#a = apply(sen,2,function(x){return(cor(x,sen[,325],use="complete.obs"))})
#t = sort(a,decreasing = TRUE)
#s = order(a,decreasing = TRUE)


## Classifier
i = 325
#i=1429
Y = sen[,i]
Y = Y[!is.na(sen[,i])]
X = GE[!is.na(sen[,i]),] # remove cell lines that are "NA" For each drug 

#Normalization
X = scale(X)
Y = scale(Y)
#Y = (Y-min(Y))/(max(Y)-min(Y))
Rep = 5
Final_Cor = rep(0,Rep)
MSE = rep(0,Rep)
for (j in 1:Rep){
  print(j)
  
  ## Split data into train & test
  sample = sample.split(Y, SplitRatio = .8)
  
  Xtrain = subset(X, sample == TRUE)
  Xtest  = subset(X, sample == FALSE)
  Ytrain = subset(Y, sample == TRUE)
  Ytest  = subset(Y, sample == FALSE)
  
  ###Model
  # Initialize a sequential model
  model = keras_model_sequential()
  
  # Add layers to the model
  model %>% 
    layer_dense(units = 100, activation = 'sigmoid', input_shape = ncol(Xtrain)) %>% 
                #kernel_regularizer = regularizer_l2(0.001)) %>%  
                #layer_dropout(rate = 0.4) %>% 
                #kernel_initializer = "normal" ,
                #bias_initializer = "Zeros") %>%  

    layer_dense(units = 10, activation = 'sigmoid') %>% 
    #layer_dropout(rate = 0.3) %>% 
    
    #layer_dense(units = 125, activation = 'relu') %>% 
    #layer_dropout(rate = 0.2) %>% 
    
    
    #layer_dense(units = 125, activation = 'relu') %>%
    layer_dense(units = 1, activation = 'linear')
  
  
  # Define an optimizer
  # sgd <- optimizer_sgd(lr = 0.01)
  
  # Compile the model
  model %>% compile(
    loss = 'mse',
    optimizer = optimizer_adam(learning_rate = 0.0001))
  
  # Fit the model 
    model %>% fit(
    Xtrain, 
    Ytrain, 
    epochs = 20, 
    batch_size = 10, 
    validation_split = 0.2)
  
  # Evaluate on test data and labels
  #score <- model %>% evaluate(Xtest, Ytest_cat, batch_size = 10)
  y_hat <- model %>% predict(Xtest)
  Final_Cor[j] = cor(Ytest,y_hat)
  MSE[j] = mean((Ytest-y_hat)^2)
  
  print(Final_Cor)
  print(MSE)
  
}
  
print(mean(Final_Cor))
#print(sd(Final_Cor))
#print(mean(MSE))
  
  
  
  