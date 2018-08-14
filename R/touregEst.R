#####  replaced by toureg.fit


touregEst <- function(x, y, start.beta, start.delta, ...)
{
    X <- x
    ini <- c(start.beta,start.delta)
    names(ini)[length(ini)] <- 'delta'
    
    out <- optim(par=ini, fn=ll, gr = gll, x=X, y=y, 
           method = 'BFGS', 
           control = list(fnscale=-1), hessian = TRUE)

    mle <- out$par  # MLE
    beta <- mle[-length(mle)]
 
    yhat <- as.vector(exp(X %*% beta))
    out$fitted.values <- yhat
  
    V <- MASS::ginv(-out$hessian) 
    SE <- sqrt(diag(V))
    names(SE) <- names(mle)     
      
    return( list(coefficients = beta,
                 delta = mle['delta'], 
                 se = SE, 
                 vcov = V, 
                 fitted.values = yhat, 
                 df = nrow(X) - ncol(X), 
                 logLik= out$value,
                 convergence = out$convergence,
                 formula = formula,
                 start.beta = start.beta, 
                 start.delta = start.delta)
     ) 
}


