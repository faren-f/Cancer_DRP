#                    Created on Wed Aug 10 17:10 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: 
library(igraph)
Infogenes = function(MyGraph,my_genes,DT){

  InitialScores = setNames(rep(0,length(my_genes)), my_genes)
  InitialScores[DT] = 1
  
  Score_diff = page_rank(MyGraph,directed = FALSE, damping = 0.99,
                         personalized = InitialScores)$vector
  Score_diff = Score_diff/max(Score_diff)


  return(Score_diff)
}



