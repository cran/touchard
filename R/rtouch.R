r0 <- function(p0, lambda, delta, N=NULL, eps=sqrt(.Machine$double.eps)){  
  # p0 = 1/tau is left out so it is not computed in every call to r0
  q <- p <- p0
  j = 0 
  u = runif(1)
  while(u > q){
     j = j + 1 
     p = p * lambda * ((j+1)^delta) / j^(delta+1)
     q = q + p
  }   
  return(j)
}

#  needed in rtouch
TAU <- function(pars) tau(lambda=pars[1], delta=pars[2], N=NULL, eps=sqrt(.Machine$double.eps) )
R0 <- function(pars) r0(p0=pars[1], lambda=pars[2], delta=pars[3], N=NULL, eps=sqrt(.Machine$double.eps))
 
rtouch <- function(n, lambda, delta, N=NULL, eps=sqrt(.Machine$double.eps)){
   if(length(lambda) == 1 & length(delta) == 1){ 
     stopifnot(lambda >= 0)
     p0 <- 1/tau(lambda,delta, N, eps)  # avoid computing this inside 'replicate'
     x <- replicate(n, r0(p0, lambda, delta, N, eps))
   }else{             # recycling as in other r-type functions 
	 stopifnot(all(lambda >= 0))
	 lam <- rep(lambda, length=n)
	 del <- rep(delta, length=n)
	 pars <- cbind(lam,del)
	 p0 <- 1/apply(pars, 1, FUN=TAU)
	 x <- apply(cbind(p0,pars), 1, FUN=R0)
   }	
   return(x)
}
