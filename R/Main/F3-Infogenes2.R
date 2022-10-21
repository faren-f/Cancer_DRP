#                    Created on Wed Aug 10 17:10 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function receives Xtrain, ytrain data, number of high correlated 
# genes (N1) with ytrain and number of high scored genes (N2) after diffusion in 
#the ppi network, to select informative genes.

Infogenes = function(Xtrain,ytrain,MyGraph,N_feat=500){
  
  Corr = abs(cor(Xtrain,ytrain))
  
  InitialScores = Corr
  Score_diff = page_rank(MyGraph,directed = FALSE,damping = 0.8,
                         personalized = InitialScores)$vector
  Score_diff = Score_diff/max(Score_diff)
  
  S_ind = order(Score_diff, decreasing = TRUE)
  good_genes = S_ind[1:N_feat]
  
  Xtrain = Xtrain[,good_genes]
  return(Xtrain)
}









