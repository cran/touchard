qnormal.toureg <- function (object, nsim = 20, level = 0.95, envelope = FALSE) 
{
    #y <- model.response(model.frame(object))
    terms <- delete.response(terms(object))
    x <- model.matrix(terms, data=object$data) 
    n <- NROW(x)
    alpha <- (1 - level)/2
    lamhat <- fitted(object)
    dhat <- object$delta 
    h <- gleverage(object)
    res <- residuals(object)/sqrt(1-h)
    if(!envelope){
		return( matrix(sort(res)) )
	}else{
    bhat <- coef(object)
    e <- matrix(0, n, nsim)
    e1 <- numeric(n)
    e2 <- numeric(n)
    Y <- matrix(
       as.numeric( rtouch(n*nsim, lambda=lamhat, delta=dhat) ),
       nrow = n ) # cada coluna eh uma amostra de tamanho n
    res.sim <- function(y){
	   suppressWarnings( 
	      foo <- toureg.fit(as.matrix(x), y, start.beta=bhat, start.delta=dhat)
	      )
  S <- diag(foo$var)
  k <- length(foo$coefficients)
  V <- foo$vcov[1:k,1:k]  
  .h <- diag( S %*% x %*% V %*% t(x) )    # hatvalues
  
  l0 <- foo$fitted.values
  d0 <- foo$delta
  
  c <- cov_yw(l0, d0)
  d <- 1 / (1+y)
  S <- diag(foo$var)

  A <- cbind(S%*%x, c)
  B <- t(cbind(x,d))
   
  invF <- foo$vcov
  
  .h <- diag( A %*% invF %*% B )    # gleverage
  
     sort(
	      residuals(foo)/sqrt(1-.h) 
	      )
  
    }
    e <- apply(Y, 2, res.sim) # cada coluna eh um vetor de n residuos
    e1 <- apply(e, 1, function(u) quantile(u,alpha))
    e2 <- apply(e, 1, function(u) quantile(u,1-alpha))
    e0 <- apply(e, 1, median)
    return( cbind(sort(res), e0, e1, e2) )
	}
}
