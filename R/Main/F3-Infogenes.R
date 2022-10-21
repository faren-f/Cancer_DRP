#                    Created on Wed Aug 10 17:10 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function receives Xtrain, ytrain data, number of high correlated 
# genes (N1) with ytrain and number of high scored genes (N2) after diffusion in 
#the ppi network, to select informative genes.

Infogenes = function(Xtrain,ytrain,MyGraph,my_genes,N_feat = 100){
  
  N2 = N_feat
  N1 = floor(N_feat*0.75)
  Corr = abs(cor(Xtrain,ytrain))
  Cor_ind = order(Corr,decreasing = TRUE)
  
  inf_gene_ind = Cor_ind[1:N1]
  inf_gene = my_genes[inf_gene_ind]
  
  InitialScores = setNames(rep(0,length(my_genes)), my_genes)
  InitialScores[inf_gene]=1
  
  Score_diff = page_rank(MyGraph,directed = FALSE,damping = 0.8,
                         personalized = InitialScores)$vector
  Score_diff = Score_diff/max(Score_diff)
  
  S_ind = order(Score_diff, decreasing = TRUE)
  S_ind = S_ind[1:N2]
  
  good_genes = my_genes[S_ind]
  
  Xtrain = Xtrain[,good_genes]
  return(Xtrain)
}









