# given mu and delta, find lambda

# using uniroot works well for lambda less than 50 and |delta| < 15 with errors around 1e-5

mdif <- function(z, m, dd, N, eps) mu(lambda=z, delta=dd, N=N, eps=eps)- m 

# foo <- function(z, m, dd, N, eps) tau(lambda=z, delta=dd, N=N, eps=eps)*(m+1) - tau(lambda=z, delta=dd+1, N=N, eps=eps)
# foo <- Vectorize(foo, "z")

lam.root <- function(MU, delta, N=100, eps=1e-6){
  stopifnot(MU >= 0)
  if(MU == 0) return(1e-8)
  if(delta == 0){ 
	   val <- MU
	  }else{
	   if(MU > 4*abs(delta)){
		     val <- MU - delta
		 }else{ 
       lo <- ifelse(delta<0, MU, max(1e-16,MU-delta))
       up <- ifelse(delta<0, MU-4*delta, min(MU, lo+0.1))
       #val <- root.sec(m=MU, d=delta, N=N, eps=eps, x0=lo, x1=up)
       #val <- root.new(m=MU, d=delta, N=N, eps=eps, x0=(up-lo)/2)
       val <- uniroot(mdif, m=MU, dd=delta, N=N, eps=eps, lower=lo, upper=up, extendInt="upX")$root   # default tol = 0.0001220703
       val <- ifelse(val < 0, 1e-8, val)
       # muito mais demorado e com erros maiores
       #LAM <- seq(lo,up,length=500)
       #DIF <-  abs( mu(LAM, d=delta, N=N, eps=eps) - MU )
       #val <- LAM[which.min(DIF)]
   }
   }
  return( val )
  }


Lambda <- Vectorize(lam.root, "MU")


