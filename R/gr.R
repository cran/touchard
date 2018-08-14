# Gradient of ll 

gll <- function(parm, x, y, N=NULL, eps=sqrt(.Machine$double.eps)){
  p <- ncol(x)
  b <- parm[1:p]
  del <- parm[p+1]
  lam <- exp(x %*% b)
  m <- mu(lam, del, N=N, eps=eps)
  gb <- t(x) %*% (y-m)
  gd <- sum( log(y+1) - kapa(lam, del, N=N, eps=eps) ) 
  return(c(gb, gd))
  }

