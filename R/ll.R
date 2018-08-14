#  loglikelihood

ll <- function(parm, x, y, N=NULL, eps=sqrt(.Machine$double.eps)){
  p <- ncol(x)
  b <- parm[1:p]
  lam <- exp(x %*% b)
  del <- parm[p+1]
  lnum <- del * log(y+1) + y * log(lam)
  lden <- lgamma(y+1) + log(tau(lam, del, N, eps))
  val <- sum( lnum-lden )
  return(val)
  }
