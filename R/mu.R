 mu0 <- function(lambda, delta, N=NULL, eps=sqrt(.Machine$double.eps)){
        if (lambda <= 0) stop("valid only for lambda > 0")
        t0        <- tau(lambda,delta,N,eps)
        m0        <- tau(lambda,delta+1,N,eps)/t0
        return(m0 - 1)
}
   #val <- sapply(lambda, mu0, delta=delta, N=N, eps=eps) 
  
mu <- Vectorize(mu0, vectorize.args = c("lambda", "delta"))
      
     
