# Gradient of ll 

 
gll2 <- function(parm, x, y, N=100, eps=sqrt(.Machine$double.eps)){
  p <- ncol(x)
  b <- parm[1:p]
  del <- parm[p+1]
  m <- c(exp(x %*% b))
  lam <- Lambda(m, del, N=N, eps=eps)
  
  s2 <- varTou(lam, del, N=N, eps=eps)
  kk <- kapa(lam, del, N=N, eps=eps)
  cv <- cov_yw(lam, del)
  
  gb <- t(x) %*% matrix((y-m)*m/s2, ncol=1)
  gd <- sum( -cv*(y-m)/s2 + log(y+1) - kk ) 
  return(c(gb,gd))
  }

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


gLL2 <- function(parm, x, y, N=100, eps=sqrt(.Machine$double.eps)){
  p <- ncol(x)
  b <- parm[1:p]
  del <- parm[p+1]
  m <- c(exp(x %*% b))
  lam <- Lambda(m, del, N=N, eps=eps)
  
  s2 <- varTou(lam, del, N=N, eps=eps)
  kk <- kapa(lam, del, N=N, eps=eps)
  cv <- cov_yw(lam, del)
  
  gb <- t(x) %*% matrix((y-m)*m/s2, ncol=1)
  gd <- sum( -cv*(y-m)/s2 + log(y+1) - kk ) 
  return(c(gb,gd))
  }


gll.beta <- function(beta, delta, x, y, N=NULL, eps=sqrt(.Machine$double.eps)){
  b <- beta
  del <- delta
  m <- c(exp(x %*% b))
  lam <- Lambda(m, del, N=N, eps=eps)
  s2 <- varTou(lam, del, N=N, eps=eps)
  gb <- t(x) %*% matrix((y-m)*m/s2, ncol=1)
  return(c(gb))
  }
