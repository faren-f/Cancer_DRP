#                  Created on Wed Sep 25 18:31 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function receives train and test data and hyper parameters 
# of the model to compute output prediction using linear regression model 
# regularized by L2 norm popular as Ridge regression model.

library(glmnet)
library(caret)

Ridge = function(ytrain, Xtrain, yval, Xval){
  
  train_data = cbind(Xtrain,ytrain)

  tune = expand.grid(alpha = 0,lambda = seq(1,20,by = 0.5))
  Result = c()
  for(i in 1:nrow(tune)){
    print(i)
    model = caret::train(ytrain ~., data = train_data,
                         method = "glmnet",
                         metric="RMSE",
                         allowParallel = TRUE,
                         tuneGrid = tune[20,])
    
    y_pred_val = predict(model,Xval)
    y_pred_test = predict(model,Xtest)
    
    Cor_val = cor(yval,y_pred_val)
    Cor_test = cor(ytest,y_pred_test)
    Pval_val = wilcox.test(y_pred_val[yval==1], y_pred_val[yval==2], alternative ="greater")$p.value
    Pval_val = -log10(Pval_val)
    Pval_test = wilcox.test(y_pred_test[ytest==1], y_pred_test[ytest==2], alternative ="greater")$p.value
    
    
    Result = rbind(Result,cbind(Cor_val,Pval_val))
}
  plot(Result[,2])
  tune[which.max(Result[,2]),2]
  return(Result)
}



