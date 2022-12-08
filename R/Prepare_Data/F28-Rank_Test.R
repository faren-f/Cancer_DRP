

rank.test = function(dist, x, alternative = "two.sided"){
  library(mnda)
  y = c(x, dist)
  
  if (alternative == "greater"){
    r = Rank(y, decreasing = TRUE)
    p.val = r[1]/length(y)
  }else if (alternative == "less"){
    r = 1/Rank(y, decreasing = FALSE)
    p.val = r[1]/length(y)
  }else if (alternative == "two.sided"){
    r1 = Rank(y, decreasing = FALSE)
    r2 = Rank(y, decreasing = TRUE)
    p.val = min(c(r1[1],r2[1])) / (length(y)/2)
  }
  
  return(p.val)
}