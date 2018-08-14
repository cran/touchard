varTou <- function(lambda, delta, N=NULL, eps=sqrt(.Machine$double.eps)){
        if (any(lambda <= 0)) 
            stop("valid only for lambda > 0")
        t0 <- tau(lambda, delta, N, eps)
        t1 <- tau(lambda, delta + 1, N, eps)/t0
        t2 <- tau(lambda, delta + 2, N, eps)/t0
        return(t2-t1^2)
}
        
#varTou <- Vectorize(v0, vectorize.args = c("lambda", "delta"))
