touqp <- function(beta, delta, x, y, N=50, eps=1e-8, dm=NULL, method = c("qp1", "qp2"), maxit){

   poi <- glm.fit( x=x, y=y, family=poisson(link='log'), start=beta, control = list(maxit=maxit) )
   m <- as.numeric( fitted(poi) )
   method <- match.arg(method)		   
 
 if(method=="qp1") {
   del.qpt <- mean(m-(y-m)^2) 
   ltil <- Lambda(MU=m, delta=del.qpt)
   v <- varTou( ltil, delta = del.qpt )
   V <- diag(v)
   M <- diag(m)
   A <- t(x) %*% M %*% x
   B <- t(x) %*% V %*% x 
   Binv <- MASS::ginv(B)
   H <- -A  %*% Binv %*% A  # for compatibility with output from optim / see toureg.fit
}

 if(method=="qp2") {
	 foo <- function(d){
       lam <- Lambda(MU=m,delta=d)
       return( sum( (y-m)^2 + (m+1)^2 - tau(lam,d+2)/tau(lam,d) ))
     }
      if(dm==0) stop("qp2 method: dm=0 implies zero-length interval to be searched for delta!")
      dm <- abs(dm)
      interv <- delta + c(-1,1)*dm
      del.qpt <- uniroot(foo, interv)$root   
      ltil <- Lambda(MU=m, delta=del.qpt)
      v <- varTou( ltil, delta = del.qpt )
      V <- diag(v)
      M <- diag(m)
      A <- t(x) %*% M %*% x
      B <- t(x) %*% V %*% x 
      Binv <- MASS::ginv(B)
      H <- -A  %*% Binv %*% A  # for compatibility with output from optim / see toureg.fit
 }
  est <- c(coef(poi),delta=del.qpt)
  f <- sum(dtouch(y, lambda=ltil, delta=del.qpt, log=TRUE))
  g <- NULL
    
  return(list(par=est, value=f, grad = g, hessian=H,  w.iwls=poi$weights, niter = poi$iter, convergence=poi$converged))
}


