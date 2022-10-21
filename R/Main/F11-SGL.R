
#                  Created on Fri Aug 19 13:55 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function receives train and test and hyper parameters of the model 
# to compute output prediction using Sparse Group Lasso(SGL)

library(SGL)
My_SGL = function(ytrain,Xtrain,Xtest,index){
  
  data = list(x=Xtrain, y = ytrain)
  cvfit = cvSGL(data, index, type = "linear", maxit = 5, nlam = 10, nfold = 5,
               lambdas = seq(.000001,.0001,.000001))
  
  i_lambda_opt = which.min(cvfit$lldiff)
  lambda_opt = cvfit$lambdas[i_lambda_opt]
  fit = SGL(data, index, type = "linear", maxit = 10, lambdas = c(lambda_opt,0))
  
  y_pred_SGL = predictSGL(fit, Xtest, 1)
  
  return(y_pred_SGL)
}
  
  
  
  
  