### E[ log(Y+1) ]

k0 <- function(lambda,delta, N=NULL, eps=sqrt(.Machine$double.eps)){
 if(is.null(N)){
    B = lambda * 2^delta * log(2)
    re <- 1
    TAU <- tau(lambda,delta, N, eps)
    K = B/TAU
    j <- 1
    while(re > eps){
      B = lambda/(j+1) * (log(j+2)/log(j+1)) * ((j+2)/(j+1))^delta * B
      K = K + B/TAU
      j=j+1
      re <- B/K
    }
  }else{
	x <- 0:N
	lnum <- delta * log(x+1) + x * log(lambda)
    lden <- lgamma(x+1) + log(tau(lambda, delta, N, eps))
    K <- sum(log(x+1)*exp(lnum-lden))
    }
return(K)
}  

kapa <-  Vectorize(k0, c('lambda','delta'))




