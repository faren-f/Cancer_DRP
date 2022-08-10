#                    Created on Wed Aug 10 12:04 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function receives Xtrain, ytrain data and number of selected 
# genes(N), to select genes that are highly correlated with drug sensitivity (ytrain).

N = 1000
high_corr = function(Xtrain,ytrain,N){
  
  Corr = abs(cor(Xtrain,ytrain))
  high_corr = order(Corr,decreasing = TRUE)
  Xtrain = Xtrain[,high_corr[1:N]]
  return(Xtrain)
}