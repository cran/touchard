
v <- function(lambda,delta,N=NULL,eps=sqrt(.Machine$double.eps)){
        if (lambda <= 0) stop("valid only for lambda > 0")
        t0        <- tau(lambda,delta,N,eps)
        m0        <- tau(lambda,delta+1,N,eps)/t0
        tau(lambda,delta+2,N,eps)/t0 - m0^2
}



qtouch0 <- function(p, lambda, delta, N=NULL,eps=sqrt(.Machine$double.eps)){
    stopifnot( lambda >= 0 )
    
    a <- function(j) tau(lambda, delta+j, N, eps)/tau(lambda,delta,N,eps)
     
    # k = cumulant, m = moment, M = central moment 
    k1 <- mu(lambda, delta,N,eps) # = m1
	k2 <- v(lambda,delta,N,eps)   # = M1
	m2 <- a(2) - 2*a(1) + 1
	m3 <- a(3) + 3*(a(1)-a(2)) - 1
	k3 <- m3 - 3*k1*m2 + 2*k1^3  # k1 = m1 and k3 = M3
    
    z <- qnorm(p)
    foo <- z + (z^2-1)*k3/6  
    cf <- round(k1 + sqrt(k2)*foo)    # order 1 cornish-fisher; initial approximation
    if( cf < 0 ) cf <- 0
    
    y <- cf
    
    if( ptouch(y, lambda, delta,N,eps) > p ){
		if( y > 1 ){
			while( ptouch(y-1, lambda, delta,N,eps) > p & y > 0 ) y <- y - 1
		} # nao tem else 
	}else{
		while( ptouch(y, lambda, delta,N,eps) < p  ) y <- y + 1
		}	
    return(y)    
}

qtouch <- Vectorize(qtouch0, 'p')

