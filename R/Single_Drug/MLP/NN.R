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
y = sen[,i]
y = y[!is.na(sen[,i])]
X = GE[!is.na(sen[,i]),] # remove cell lines that are "NA" For each drug 

#Normalization
X = scale(X)
y = scale(y)
#y = (y-min(y))/(max(y)-min(y))
Rep = 5
corr = rep(0,Rep)
mse = rep(0,Rep)
for (j in 1:Rep){
  print(j)
  
  ## Split data into train & test
  sample = sample.split(y, SplitRatio = .8)
  
  Xtrain = subset(X, sample == TRUE)
  Xtest  = subset(X, sample == FALSE)
  ytrain = subset(y, sample == TRUE)
  ytest  = subset(y, sample == FALSE)
  
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
    ytrain, 
    epochs = 20, 
    batch_size = 10, 
    validation_split = 0.2)
  
  # Evaluate on test data and labels
  #score <- model %>% evaluate(Xtest, ytest_cat, batch_size = 10)
  y_pred <- model %>% predict(Xtest)
  corr[j] = cor(ytest,y_pred)
  mse[j] = mean((ytest-y_pred)^2)
  
  print(corr)
  print(mse)
  
}
  
print(mean(corr))
#print(sd(corr))
#print(mean(mse))
  
  
  
  