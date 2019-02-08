touch_iwls <- function(beta, delta, x, y, tol=1e-6, tmax=25, N=NULL, eps=1e-6){
 b <- beta
 gb <- gll.beta(b, delta=delta, x=x, y=y)
 t <- 1
 while( max(abs(gb)) > tol & t <= tmax){
  eta <- drop( x %*% b )
  m <- exp( eta )
  lam <- Lambda(MU=m, delta=delta, N=N, eps=eps)
  s2 <- varTou(lambda=lam, delta=delta, N=N, eps=eps)
  w <- drop( m^2/s2 )
  q <- m*(y-m)/s2
  invfish <- MASS::ginv( crossprod(x, w*x) )
  b <- b + tcrossprod(invfish,x) %*% q
  gb <- gll.beta(b, delta=delta, x=x, y=y)
  t <- t+1
 }
 maxl <- ll2(beta=b, delta=delta, x=x, y=y, N=N, eps=eps)
 return(list(beta = b, grad = gb, value=maxl, w=w, niter = t-1))
}
