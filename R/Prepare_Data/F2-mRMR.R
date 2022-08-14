
#                  Created on Thu Jul 28 16:16 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function find features that have highest corelation with drug sensitivity 
# while they have lowest correlation with other previously selected features.



mRMR = function(Xtrain, ytrain, N_feat = NA, alpha=1, do.plot = TRUE){
  
  if (is.na(N_feat)){
    N_feat = ncol(Xtrain)
  
  }else{
    cor_xy = abs(cor(Xtrain,ytrain))
    Cor_xx = abs(cor(Xtrain,Xtrain))
    
    cor_xy_orders = order(cor_xy, decreasing = TRUE)
    selected_feat = cor_xy_orders[1]
    residual_feat = cor_xy_orders[-1]
    scores = max(cor_xy)
    
    while(length(residual_feat)!= (ncol(Xtrain)-N_feat)){
      S = c()
      for(j in residual_feat){
        score_j = cor_xy[j] - alpha*(sum(Cor_xx[j, selected_feat])/length(selected_feat))
        S = c(S, score_j)
      }
      selected_feat = c(selected_feat, residual_feat[which.max(S)])
      residual_feat = residual_feat[-which.max(S)]
      scores = c(scores, max(S))
    }
  }
  Xtrain = Xtrain[,selected_feat]
  if (do.plot)
    plot(scores)
  return(Xtrain)
  }



