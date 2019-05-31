residuals.toureg <- function (object, type = c("deviance",  "pearson", "response"), N = 100, eps = 1e-6, ...) 
{
    res <- object$residuals
    type <- match.arg(type)
    
    if (type == "response")  return(res)
    
    if (type == "pearson")  return( res/sqrt(object$var) )
   
    if (type == "deviance"){
     if(object$method == "qp1" || object$method=="qp2")
        return(NULL) 
        else 
         y <- if (is.null(object$y)){
	      model.response(model.frame(object)) 
	         }else{
		    object$y
		}
     	   ltilde <- Lambda(MU=y, delta=object$delta, N=N, eps=eps)
           lhat <- object$lambda
	   tau.tilde <- tau(lambda=ltilde, delta=object$delta, N=N, eps=eps)
	   tau.hat <- tau(lambda=lhat, delta=object$delta, N=N, eps=eps)
	       return( sign(res)*sqrt( 2*(y*(log(ltilde) - log(lhat)) - log(tau.tilde) + log(tau.hat)) ) )
     }
}
  
