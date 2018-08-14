score.delta <- function(x, freq=NULL, max=50, data){

  max1 <- max+1
# Score test
 if(class(x) == "formula"){
  data <- na.omit(data)
  mf <- model.frame(formula=x, data=data)
  terms <- attr(mf, "terms")
  X = model.matrix(terms, mf)
  y = model.response(mf)
  z <- log(y+1)
  Sz <- sum(z)
  ltil <- fitted( glm(x, data=data, family=poisson) )
  Sk <- sum( kapa(ltil,0) )

  vz <- function(lam) {
	mz <- sum(log(1:max1) * dpois(0:max,lam))
	sum( (log(1:max1)-mz)^2 * dpois(0:max,lam) )
  }
  
  Sv <- sum( sapply(ltil, vz) )

  cv <- function(lam) {
	y <- 0:max
	z <- log(y+1)
	p <- dpois(y,lam)
	ey <- lam
	ez <- sum(z*p)
	sum((y-ey)*(z-ez)*p)
	}
	
  g <-  sapply(ltil, cv) 
  dim(g) <- c(length(y),1)
  x.g <- crossprod(X,g)
  inv <- solve( crossprod(X, ltil * X) )
  cx <- t(x.g) %*% crossprod(inv, x.g)   # inv is symmetric so inv'xg = inv gx
   
  S <- (Sz - Sk)^2 / (Sv - cx)
 
}else{

 if(is.table(x)){
	 if(!is.null(freq)) warning("_freq_ has been ignored since _x_ is already of class 'table'")
	 ny <- as.vector(x)
	 y <- rep( as.numeric(names(x)), ny )
 }else{
	if( is.null(freq) ) y <- x else y <- rep(x, freq) 
  }
  z <- log(y+1)
  n <- length(y)
  l <- mean(y)
  k <- kapa(l,0)
  p <- dpois(0:max,l)
  mz <- sum(log(1:max1) * p)
  v <- sum( (log(1:max1)-mz)^2 * p )
  ez <- sum(log(1:max1)*p)
  g <- sum((0:max -l)*(log(1:max1)-ez)*p)
  s2 <- varTou(l,0)  # = l

  num <- n*(mean(z) - k)^2
  den <- v - g^2/s2
  S <- num/den
  
 }
 
 return(list(stat=S, pval=1-pchisq(S, df=1)))
}
