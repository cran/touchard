predict.toureg <- function(object, newdata=NULL, type=c("response", "lambda", "variance", "linpred"), se.fit=NULL, ...){	
	
	 type <- match.arg(type)
     terms <- delete.response(terms(object))
     
     if(is.null(newdata)){
	  lamhat <- object$lambda
      muhat <- object$mu
      v <- object$var
	  Z <- model.matrix(terms, data=object$data) 
	  linpred <- as.vector(Z %*% coef(object))
        }else{
	# formula interface
	  Z <- model.matrix(terms, model.frame(terms, as.data.frame(newdata)))
	  linpred <- as.vector(Z %*% coef(object))
	  lamhat <- switch(object$regress, 
	                   lambda = exp( linpred ),
	                   mu = Lambda( exp( linpred ), object$delta, N=NULL )
	                   )
	  muhat <- switch(object$regress, 
	                   lambda = mu( exp( linpred ), object$delta, N=NULL ),
 	                   mu = exp( linpred )
	                   )
	  v <- varTou( lamhat, object$delta, N=NULL )
	                   
	}

    p <- ncol(Z)	
    S <- object$vcov[1:p,1:p]   # vcov is Var[ (b1,..,bp,delta) ]
    n <- nrow(Z)
     
    if(object$regress == "lambda"){
    if(is.null(se.fit)) { ep <- NULL }else{ 
	  ep <- switch(type, 	
	               lambda = lamhat*sqrt( diag( Z %*% S %*% t(Z) ) ),  
	               response = v*sqrt( diag( Z %*% S %*% t(Z) ) ),    
	               variance = NULL,
	               linpred = NULL) 
	  }
    
      out <- switch(type, 
                    lambda = lamhat,
                    response = muhat,
                    variance = v,
                    linpred = linpred)
	
    
}
	if(object$regress == "mu"){
    if(is.null(se.fit)) { ep <- NULL }else{ 
	  ep <- switch(type, 	
	               lambda = (lamhat*muhat/v) * sqrt( diag( Z %*% S %*% t(Z) ) ),  
	               response = muhat*sqrt( diag( Z %*% S %*% t(Z) ) ),    
	               variance = NULL,
	               linpred = NULL) 
	  }
    
      out <- switch(type, 
                    lambda = lamhat,
                    response = muhat,
                    variance = v,
                    linpred = linpred)
	
}
	 return(list(fit=out, se.fit=ep))

}

