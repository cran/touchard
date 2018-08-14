
hatvalues.toureg <- function(object, ...){
  sig2 <- object$var
  S <- diag(sig2)
  k <- length(object$coefficients)
  V <- object$vcov[1:k,1:k]
  
  terms <- delete.response(terms(object))
  x <- if (is.null(object$x)){
	   model.matrix(terms, model.frame(object))
	   }else{
		    object$x
		}
    
  h <- S %*% x %*% V %*% t(x) 
  return(diag(h))
}
