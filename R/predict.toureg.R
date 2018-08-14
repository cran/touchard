predict.toureg <- function(object, newdata=NULL, type=c("invlink", "response", "variance"), se.fit=NULL, ...){	
	
	 type <- match.arg(type)
     terms <- delete.response(terms(object))
     
     if(is.null(newdata)){
	  lamhat <- fitted(object)
	  Z <- model.frame(terms, data=object$data)
        }else{
	# formula interface
	  Z <- model.matrix(terms, model.frame(terms, as.data.frame(newdata)))
	  lamhat <- exp( as.vector(Z %*% coef(object)) )
	}

    p <- ncol(Z)	
    S <- object$vcov[1:p,1:p]   # vcov is Var[ (b1,..,bp,delta) ]
    n <- nrow(Z)
    v <- varTou(lamhat, object$delta)

    if(is.null(se.fit)) { ep <- NULL }else{ 
	  ep <- switch(type, 	
	               invlink = lamhat*sqrt( diag( Z %*% S %*% t(Z) ) ),  
	               response = v*sqrt( diag( Z %*% S %*% t(Z) ) ),    
	               variance = NULL) 
	  }
    
    out <- switch(type, 
                  invlink = lamhat,
                  response = mu(lamhat, object$delta),
                  variance = v)
	
    return(list(fit=out, se.fit=ep))	
}

