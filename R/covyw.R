### E[ Y x W ] and COV(Y,W), W=log(Y+1)

cov_yw <- function(lambda,delta,N=100){
 c0 <- function(lambda, delta, N){
    B = 0
    j=1
    B[j] = lambda * 2^delta * log(2)
    eps = sqrt(.Machine$double.eps)
    while(B[j] > eps){
      B[j+1] = (lambda/j) * (log(j+2)/log(j+1)) * ((j+2)/(j+1))^delta * B[j]
      j=j+1
    }
    return(sum(B)/tau(lambda,delta,N))
    }
 EYW <- sapply(lambda, c0, delta=delta, N=N) 
 EY <- mu(lambda, delta, N=N)
 EW <- kapa(lambda, delta, N=N)
 val <- EYW - EY*EW
 return(val)
}  

