#                  Created on Wed Dec 8 15:07 2022

#                     @author: Farzaneh Firoozbakht
#                   Email: faren.firoozbakht@gmail.com

# Description: This function receives a value of Beta and a vector of Beta_Null
# which is the Null distribution, and return the probability value of the 
# extremeness of Beta with diffrent alternatives of "less", "greater" and "two.sided"


Probability_norm = function(Beta, Beta_Null, alternative){
  
  if(alternative == "less"){
    p.val = pnorm(Beta, mean(Beta_Null), sd(Beta_Null), lower.tail = TRUE)
  }else if(alternative == "greater"){
    p.val = pnorm(Beta, mean(Beta_Null), sd(Beta_Null), lower.tail = FALSE)
  }else if(alternative == "two.sided"){
    p.val = 2 * pnorm(abs(Beta), mean(Beta_Null), sd(Beta_Null), lower.tail = FALSE)
  }
  return(p.val)
}



