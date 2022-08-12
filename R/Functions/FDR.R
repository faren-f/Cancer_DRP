FDR = function(x,y){
  
  # Fisher discriminant ratio(FDR) calculation.
  # Reference: Theodoridise, machine learning 
  # Author: Farzaneh Firoozbakht
  # x is a feature matrix [sample x dimention]
  # y is a vector of labels 
  
  assertthat::assert_that(nrow(x)==length(y))
  
  x = data.frame(x)
  class = unique(y)
  
  P_classes = c()
  for(i in class){
    No_classes = sum(y == i)
    P_classes = c(P_classes, No_classes/length(y))
  }
  
  ## Within class scatter matrix
  Sw = 0
  for(i in 1:length(class)){
    xi = x[y==class[i],]
    Si = ifelse(ncol(x)==1, var(xi), cov(xi))
    Sw = Sw + (P_classes[i]* Si)
  }

  ## Mixture scatter matrix
  Sm = ifelse(ncol(x)==1, var(x), cov(x))
  ## FDR
  J = (sum(diag(Sm)))/(sum(diag(Sw)))
  
  return(J)
}

