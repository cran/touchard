touglm <- function(beta, delta, x, y, 
                       etol=sqrt(.Machine$double.eps), gtol=1e-4, reltol=sqrt(.Machine$double.eps), 
                       maxit=20, N=50, eps=1e-8, dm=10){
						   
 if(dm==0) stop("glm method: dm=0 implies zero-length interval to be searched for delta!")
 dm <- abs(dm)
 b <- beta
 d <- delta
 f <- ll2(beta=b, delta=d, x=x, y=y, N=N, eps=eps)
 g <- gll2(c(b, d), x=x, y=y, N=N, eps=eps)
 t <- 1
 dif.e <- 1
 dif.f <- 1
 rhs.e <- 0
 rhs.f <- 0 
 
 while( dif.e > rhs.e & max(abs(g)) > gtol & dif.f > rhs.f & t <= maxit ){
	old.e <- c(b,d)
	old.f <- ll2(beta=b, delta=d, x=x, y=y, N=N, eps=eps)
	
    out1 <- touch_iwls(beta=b, delta=d, x=x, y=y, tol=gtol, tmax=maxit, N=N, eps=eps) 
    if(out1$niter==maxit) warning("IWLS step reached _maxit_")
    b <- out1$beta	
    interv <- d + c(-1,1)*dm
    out2 <- optimize(ll2, interval=interv, maximum=TRUE, beta=b, x=x, y=y, N=N, eps=eps)
    d <- out2$maximum
    g <- gll2(c(b,d), x=x, y=y, N=N, eps=eps)
    t <- t+1
    
    dif.e <- crossprod(c(b,d)-old.e)
    rhs.e <- etol*(etol+crossprod(c(b,d)))
    f <- ll2(beta=b, delta=d, x=x, y=y, N=N, eps=eps)
    dif.f <- abs(f-old.f)
    rhs.f <- reltol*(reltol+abs(f))    
    
}
  X <- x
  est <- c(b,d)
  FUN <- function(x) LL2(parm=x, x=X, y=y, N=N, eps=eps)
  H <- numDeriv::hessian(func=FUN, x=est, method="Richardson")
  
  conv <- ifelse(t==1+maxit, 1L, 0L)
  
  return(list(par=est, value=f, grad = g, hessian=H,  w.iwls=out1$w, niter = t-1, convergence=conv))
}


