#lf <- function(parms, x) sum( dtouch(x, lambda=parms[1], delta=parms[2], log=TRUE) )
lf <- function(parms, x, freq = NULL, rc=FALSE, trunc.at.zero=FALSE){
	lam <- parms[1]
	d <- parms[2]
	if(is.null(freq)){ 
		val <- sum( dtouch(x, lambda=lam, delta=d, log=TRUE) ) 
		}else{
		if(length(x) != length(freq)) stop("x and freq must have the same length")
		g <- length(x)
		k <- g-1
		if(k < 2) stop("largest observed count should be at least 2")
		if(rc){
		  .f <- freq[-g] 
		  fg <- freq[g]
		  .x <- x[-g]
		  xg <- x[g]
		  l1 <- sum( .f*dtouch(.x, lambda=lam, delta=d, log=TRUE) )
          l2 <- fg*log((1-ptouch(xg-1,lambda=lam, delta=d)))
          val <- l1+l2
		}else{
		  val <- sum( freq*dtouch(x, lambda=lam, delta=d, log=TRUE) ) 
        }
	    if(trunc.at.zero){ 
		   l0 <- sum(freq) * log(1-dtouch(0, lambda=lam, delta=d))
		   val <- val - l0 
		}
	}
	return(val)
}

ml <- function(x, freq, start, rc, tz){
 
 if(is.null(freq)){ n <- length(x) }else{ n <- sum(freq) }
 
 fm <- optim(c(start$lambda, start$delta), lf, method='L-BFGS-B', lower = c(1e-8,-Inf), 
              hessian=TRUE, x=x, freq=freq, rc=rc, trunc.at.zero=tz, 
              control=list(fnscale=-1))
 ll <- fm$value
 AIC  <-  -2*ll+2*2
 BIC  <-  -2*ll+log(n)*2
 V <- -MASS::ginv(fm$hessian)

#   LR test
 if(is.null(freq)){ 
	 xbar <- mean(x)
	 lpois <- sum(dpois(x, lambda = xbar, log = TRUE))
 }else{
	 xbar <- weighted.mean(x, w=freq)
	 lpois <- sum(freq*dpois(x, lambda = xbar, log = TRUE))
 }
 if(tz){ lpois <- lpois - sum(freq) * log(1-dpois(0, lambda=xbar)) }
 lr.stat <- -2*(lpois-ll)   

# Wald test
 del <- fm$par[2]
 w.stat <- del^2 / V[2,2]

 c1 <- pchisq(lr.stat, df=1, lower.tail=FALSE)
 c2 <- pchisq(w.stat, df=1, lower.tail=FALSE)
 
 ep <- sqrt(diag(V))
 est <- fm$par
 names(est) <- names(ep) <- c("lambda", "delta")
 
 fit <- list(estimate=est, se=ep, vcov=V, message=fm$message) 
 
 r <- data.frame( test.stat = c(lr.stat,w.stat), pval = c(c1,c2),  row.names=c("LR","Wald") )
 
 val <- list(fit=fit, aic=AIC, bic=BIC, test=r, method="ml", data=list(x=x,freq=freq))
 return(val)
}

mm1 <- function(x, start, freq){
 
 if(is.null(freq)){ y <- x }else{ y <- rep(x,freq) }
 
 ybar <- mean(y)
 y2bar <- mean(y^2)
 
 H <- function(parms){
	lam <- parms[1]
	del <- parms[2]
	a <- function(j) tau(lam, del+j)/tau(lam,del)
	m1 <- a(1) - 1
	m2 <- a(2) - 2*a(1) + 1
	h1 <- m1 - ybar
	h2 <- m2 - y2bar
	
    return(c(h1,h2))
}

# require('nleqslv')
 MM <- nleqslv::nleqslv(c(start$lambda,start$delta), H, jacobian = TRUE)
 est <- MM$x
 
# require('numDeriv')
 G <- numDeriv::jacobian(H, est)
 Ginv <- MASS::ginv(G)

# Ginv <- solve( mm3$jac )  # nao bate com Ginv acima mas produz S abaixo similar !?

 V <- cov(cbind(y, y^2))
 
 S <- Ginv %*% V %*% t(Ginv) / length(y)
 ep <- sqrt(diag(S))
 
 names(est) <- names(ep) <- c("lambda", "delta")
 
 fm <- list(estimate = est, se = ep, vcov = S)
  
 AIC <- NULL
 BIC <- NULL
 
 # Wald test
 del <- as.numeric( est['delta'] )
 w.stat <- del^2 / S[2,2]

 
 pv <- pchisq(w.stat, df=1, lower.tail=FALSE)

 r <- data.frame( test.stat = w.stat, pval = pv, row.names="Wald"  )
 
 val <- list(fit=fm, aic=AIC, bic=BIC, test=r, method="mm", data=list(x=x,freq=freq))
 return(val)
}

mm2 <- function(x, start, freq){

 if(is.null(freq)){ y <- x }else{ y <- rep(x,freq) }
 
 ybar <- mean(y)
 w <- log(1+y)
 wbar <- mean(w)
 
 H <- function(parms){
	lam <- parms[1]
	del <- parms[2]
	a <- function(j) tau(lam, del+j)/tau(lam,del)
	m1 <- a(1) - 1
	h1 <- m1 - ybar
	h2 <- kapa(lam,del) - wbar
	
    return(c(h1,h2))
}

# require('nleqslv')
 MM <- nleqslv::nleqslv(c(start$lambda,start$delta), H, jacobian = TRUE)
 est <- MM$x
 
# require('numDeriv')
 G <- numDeriv::jacobian(H, est)
 Ginv <- MASS::ginv(G)

# Ginv <- solve( mm3$jac )  #nao bate com Ginv acima mas produz S abaixo similar !?

 V <- cov(cbind(y, w))
 
 S <- Ginv %*% V %*% t(Ginv) / length(y)
 ep <- sqrt(diag(S))
 
 names(est) <- names(ep) <- c("lambda", "delta")
 
 fm <- list(estimate = est, se = ep, vcov = S)
  
 AIC <- NULL
 BIC <- NULL
 
 # Wald test
 del <- as.numeric( est['delta'] )
 w.stat <- del^2 / S[2,2]

 
 pv <- pchisq(w.stat, df=1, lower.tail=FALSE)

 r <- data.frame( test.stat = w.stat, pval = pv, row.names="Wald" )

 val <- list(fit=fm, aic=AIC, bic=BIC, test=r, method="gmm", data=list(x=x,freq=freq))
 return(val)
}


toufit <- function(x, freq = NULL, start, method = c("ml","mm","gmm"), rc=FALSE, trunc.at.zero=FALSE){
 
 if(is.table(x)){
	 if(!is.null(freq)) warning("_freq_ has been ignored since _x_ is already of class 'table'")
	 ny <- as.vector(x)
	 y <- as.numeric(names(x))
 }else{
	 if(is.null(freq)){ 
		   X <- table(x)
		   ny <- as.vector(X)
	       y <- as.numeric(names(X))
		}else{ 
		   ny <- freq  
		   y <- x
			 } 
 
	 }
	 
 method <- match.arg(method)
 if( any(y < 0) ) stop("counts must be nonnegative")
 if( any(ny < 0) ) stop("frequencies must be nonnegative")
 if( method != "ml" & any(c(rc,trunc.at.zero)) ) stop("Method of Moments not available for either censored or truncated data.")
 
 if( trunc.at.zero & any(x==0) ) stop("truncation at zero set to 'TRUE' but data contain zeroes")
  
 if( missing(start) ){
    
    n <- sum(ny)
    lhs <- lgamma(y+1) + log(ny) - log(n)
    ly1 <- log(y+1) 
    myline <- coef(lm(lhs~y+ly1))
    l0 <- exp(myline[2])
    d0 <- myline[3]
       
    start <- list( lambda=l0, delta=d0 ) 
 }else{
	 start <- as.list(start)
	 if(is.null(names(start))) stop("_start_ must be list(lambda=..., delta=...)")
	 parnames <- sort(names(start))
	 if( !all(parnames == c("delta", "lambda")) ) stop("_start_ must be list(lambda=..., delta=...)")
	 }

 val <- switch(method, 
			   ml = ml(x=y, start=start, freq=ny, rc=rc, tz=trunc.at.zero),
			   mm = mm1(x=y, start=start, freq=ny),
			   gmm = mm2(x=y, start=start, freq=ny)
			   )
 class(val) <- "toufit" 
 val$trunc.at.zero <- trunc.at.zero
 val$rc <- rc
 val$start <- start
 
 return( val )
 }



