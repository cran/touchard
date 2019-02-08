tau0 <- function(lambda,delta,N=NULL,eps=sqrt(.Machine$double.eps)) {
  a <- lambda
  b <- delta
  if(a < 0) stop("tau cannot compute with negative lambda")
  #stopstopifnot(a >= 0)
  
  if( is.null(N) ){
     A <- 1
     TAU <- 1
     j <- 0
     re <- 1
     while(re > eps){
       k <- (a/(j+1))*((j+2)/(j+1))^b
       A <- A*k
       TAU <- TAU + A
       j <- j+1
       re <- A/TAU 
     }
  }else{
     j <- 0:N
     lnum <- b * log(j+1) + j * log(a)
     lden <- lgamma(j+1)
     TAU <- sum( exp(lnum-lden) )
  }
  return( TAU )
}

tau <- Vectorize(tau0, c('lambda','delta'))
