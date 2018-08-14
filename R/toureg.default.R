toureg.default <- function(x, y, start.beta, start.delta, ...)
{
  if( any(y<0) ) stop("response must be (nonnegative) counts")
  
  if( missing(start.beta) ){
	b0 <- coef( glm.fit(x, y, family = poisson(link = "log") ) ) 
  }else{
	b0 <- start.beta
  }

  if( missing(start.delta) ){
      d0 <- sign( mean(y) - var(y) )  # overdisp <=> delta < 0
  }else{
      d0 <- start.delta
  }


  fit <- toureg.fit(x=as.matrix(x), y=y, start.beta=b0, start.delta=d0)  
  
  fit$call <- match.call()
  return(fit)  
}


