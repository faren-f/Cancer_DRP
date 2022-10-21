#                  Created on Wed Aug 4 12:04 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function receives train and test data and hyper parameters 
# of the model to compute output prediction using MLP neural network model.

library(keras)
library(tensorflow)

MLP = function(ytrain,Xtrain,Xtest){
  model = keras_model_sequential()
  
  # Add layers to the model
  model %>% 
    layer_dense(units = 150, activation = 'sigmoid', input_shape = ncol(Xtrain),
    kernel_regularizer = regularizer_l2(0.5)) %>% 
    #layer_dropout(rate = 0.2)) %>% 
    #kernel_initializer = "normal" ,
    #bias_initializer = "Zeros") %>%  
    
    layer_dense(units = 80, activation = 'sigmoid')%>%  
    #kernel_regularizer = regularizer_l2(0.1)) %>% 
    #layer_dropout(rate = 0.3)) %>% 
    
    #layer_dense(units = 10, activation = 'sigmoid') %>% 
    #kernel_regularizer = regularizer_l2(1)) %>% 
    
    #layer_dropout(rate = 0.2)) %>% 
    
    
    #layer_dense(units = 125, activation = 'relu') %>%
    layer_dense(units = 1, activation = 'linear')
  
  
  # Define an optimizer
  # sgd <- optimizer_sgd(lr = 0.01)
  
  # Compile the model
  
  model %>% compile(
    loss = 'mse',
    optimizer = optimizer_adam(learning_rate = 0.0001))
  
  
  callbacks = list(callback_early_stopping(monitor = "val_loss", patience = 5, 
                                           restore_best_weights = TRUE))
  
  # Fit the model 
  model %>% fit(
    Xtrain, 
    ytrain, 
    epochs = 200, 
    batch_size = 5, 
    validation_split = 0.2,
    callbacks = callbacks)
  
  y_pred = predict(model, Xtest)
  
  return(y_pred)
  
}
