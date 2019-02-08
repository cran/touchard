
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
    	
  h <- hvalues(object)
  
  lhat <- object$lambda
  dhat <- object$delta
  gam <- cov_yw(lhat, dhat)
  s2 <- object$var
   
  if(object$regress == "lambda"){
	   
   d <- 1 / (1+y)
   
   S <- diag(s2)
   A <- cbind(S%*%x, gam)
   B <- t(cbind(x,d))
   
   invF <- object$vcov
  
   glev <- diag( A %*% invF %*% B )
  
  }else{  # regress == 'mu'
   
   d <-  1 / (1+y) - gam/s2
   
   D <- diag(object$mu)

   A <- cbind(D%*%x, 0)
   B <- rbind(crossprod(x, D/s2),t(d))
   
   invF <- object$vcov
  
   glev <- diag( A %*% invF %*% B )
}
   return(glev)
}
