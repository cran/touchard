toureg.fit <- function(x, y, start.beta, start.delta, 
        parscale = rep.int(1, length(start.beta)+1),  maxit = 100L, abstol = -Inf, 
        reltol = sqrt(.Machine$double.eps), N=100, eps=1e-8, ...)
{
    y <- as.vector(as.matrix(y))
    ini <- c(start.beta,start.delta)
    names(ini)[length(ini)] <- 'delta'
    
    out <- optim(par=ini, fn=ll, gr = gll, x=x, y=y, 
           method = 'BFGS', 
           control = list(fnscale=-1, parscale=parscale, maxit=maxit, abstol=abstol, reltol=reltol), 
           hessian = TRUE, N=N, eps=eps)

    mle <- out$par  # MLE
    beta <- mle[-length(mle)]
    dhat <- mle['delta']
    
    lamhat <- as.vector(exp(x %*% beta))
    out$fitted.values <- lamhat
  
    V <- MASS::ginv(-out$hessian) 
    SE <- sqrt(diag(V))
    names(SE) <- names(mle)     
    
    yhat <- mu(lambda=lamhat, delta=dhat)
    sig2 <- varTou(lambda=lamhat, delta=dhat) 
     
      
    fit <- list(coefficients = beta,
                 delta = mle['delta'], 
                 se = SE, 
                 vcov = V, 
                 fitted.values = lamhat, 
                 mu = yhat,
                 var = sig2, 
                 residuals = (y - yhat)/sqrt(sig2),  # pearson   
                 df = nrow(x) - ncol(x), 
                 logLik= out$value,
                 convergence = out$convergence,
                 start.beta = start.beta, 
                 start.delta = start.delta,
                 data=data)
                 
    class(fit) <- "toureg"
    return(fit)
}


