dtouch0 <- function (x, lambda, delta, N=NULL, eps=sqrt(.Machine$double.eps), log = FALSE){
  stopifnot(lambda >= 0)
  if( abs(x - round(x))  > .Machine$double.eps^0.5 ){ 
    warning("non-integer x = ", x)
    return(0)
  }else{ 
    lnum <- delta * log(x+1) + x * log(lambda)
    lden <- lgamma(x+1) + log(tau(lambda, delta, N, eps))
    if(log) return(lnum-lden) else return(exp(lnum-lden))
  }
}

dtouch <- Vectorize( dtouch0, c('x','lambda') )
