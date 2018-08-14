
gleverage.toureg <- function(object, ...){
  y <- if (is.null(object$y)){
	   model.response(model.frame(object)) 
	   }else{
		    object$y
		}
  terms <- delete.response(terms(object))	
  x <- if (is.null(object$x)){
	   model.matrix(terms, model.frame(object))
	   }else{
		    object$x
		}
    	
  h <- hatvalues(object)
  
  lhat <- object$fitted.values
  dhat <- object$delta
  
  c <- cov_yw(lhat, dhat)
  d <- 1 / (1+y)
  S <- diag(object$var)

  A <- cbind(S%*%x, c)
  B <- t(cbind(x,d))
   
  invF <- object$vcov
  
  glev <- diag( A %*% invF %*% B )
  return(glev)
}
