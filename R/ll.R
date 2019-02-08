#  loglikelihood

# lambda predicted
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

# mu predicted 
LL2 <- function(parm, x, y, N=NULL, eps=sqrt(.Machine$double.eps)){
  p <- ncol(x)
  b <- parm[1:p]
  del <- parm[p+1]
  .mu <- exp(x %*% b)
  lam <- Lambda(MU=.mu, delta=del, N=N, eps=eps)
  lnum <- del * log(y+1) + y * log(lam)
  lden <- lgamma(y+1) + log(tau(lam, del, N, eps))
  val <- sum( lnum-lden )
  return(val)
}
 
# versao com separacao de beta e delta  
ll2 <- function(beta, delta, x, y, N=NULL, eps=sqrt(.Machine$double.eps)){
  b <- beta
  .mu <- exp(x %*% b)
  del <- delta
  lam <- Lambda(MU=.mu, delta=del, N=N, eps=eps)
  lnum <- del * log(y+1) + y * log(lam)
  lden <- lgamma(y+1) + log(tau(lam, del, N, eps))
  val <- sum( lnum-lden )
  return(val)
  }



 
