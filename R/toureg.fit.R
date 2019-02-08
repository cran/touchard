toureg.fit <- function(x, y, start.beta, start.delta, 
        parscale = rep.int(1, length(start.beta)+1),  maxit, abstol = -Inf, 
        reltol = 1e-6, etol = 1e-6, gtol = 1e-4,
        N=100, eps=1e-4, dm = 10, regress = c("lambda", "mu"), 
        method = c( "BFGS", "CG", "Nelder-Mead", "glm", "qp1", "qp2"), ...)
{
	method <- match.arg(method)
    regress <- match.arg(regress)
    
    if(method == "Nelder-Mead") method <- "Nelder"
    
    if( missing(maxit) ){  
	
		maxit <- switch(method, 
		              glm = 10L,
		              qp1 = 25L,
		              qp2 = 25L,
                      CG = 100L,
                      BFGS = 100L,
                      Nelder =  200L
                        )
           }
  
    y <- as.vector(as.matrix(y))
    ini <- c(start.beta,start.delta)
    names(ini)[length(ini)] <- 'delta'
    
    if(regress == "mu"){
	
		out <- switch(method, 
		              glm = touglm(beta=start.beta, delta=start.delta, x=x, y=y, 
                        etol=etol, gtol=gtol, reltol=reltol, 
                        maxit=maxit, N=N, eps=eps, dm=dm),
                      CG = optim(par=ini, fn=LL2, gr=gLL2, x=x, y=y, 
                        method = 'CG', 
                        control = list(fnscale=-1, parscale=parscale, maxit=maxit, abstol=abstol, reltol=reltol), 
                        hessian = TRUE, N=N, eps=eps),
                      BFGS = optim(par=ini, fn=LL2, gr=gLL2, x=x, y=y, 
                        method = 'BFGS', 
                        control = list(fnscale=-1, parscale=parscale, maxit=maxit, abstol=abstol, reltol=reltol), 
                        hessian = TRUE, N=N, eps=eps),
                      Nelder =  optim(par=ini, fn=LL2, gr = NULL,  x=x, y=y, 
                        method = 'Nelder-Mead', 
                        control = list(fnscale=-1, parscale=parscale, maxit=maxit, abstol=abstol, reltol=reltol), 
                        hessian = TRUE, N=N, eps=eps),
                      qp1 = touqp(beta=start.beta, delta=start.delta, x=x, y=y, 
                        maxit=maxit, N=N, eps=eps, method = "qp1"),
                      qp2 = touqp(beta=start.beta, delta=start.delta, x=x, y=y, 
                        maxit=maxit, N=N, eps=eps, dm=dm, method = "qp2")
                        )
           
		mle <- out$par  # MLE
        beta <- mle[-length(mle)]
        dhat <- mle[length(mle)]
        names(mle) <- c(colnames(x), "delta")
		
		fitval <- muhat <- as.vector( exp(x %*% beta) ) 
		
		V <- MASS::ginv(-out$hessian) 
        SE <- sqrt(diag(V))
        if(method=="qp1" || method=="qp2") 
          names(SE) <- names(mle)[length(beta)]       # currently no SE for delta under QPT
        else
          names(SE) <- names(mle)     
   
    
        lamhat <- Lambda(MU=muhat, delta=dhat, N=N, eps=eps)
        sig2 <- varTou(lambda=lamhat, delta=dhat, N=N, eps=eps) 
        w <- muhat^2/sig2  # same as out$w.iwls if method == glm
    
		
       }else{  # regress="lambda"
		
		if(method == "glm" || method == "qp1" || method=="qp2") stop("not a valid method when regressing _lambda_")
		
		out <- optim(par=ini, fn=ll, gr = gll, x=x, y=y, 
           method = method, 
           control = list(fnscale=-1, parscale=parscale, maxit=maxit, abstol=abstol, reltol=reltol), 
           hessian = TRUE, N=N, eps=eps)
           
    mle <- out$par  # MLE
    beta <- mle[-length(mle)]
    dhat <- mle[length(mle)]
    
    fitval <- lamhat <- as.vector(exp(x %*% beta))
  
    V <- MASS::ginv(-out$hessian) 
    SE <- sqrt(diag(V))
    if(method=="qp") 
       names(SE) <- names(mle)[length(beta)]       # currently no SE for delta under QPT
    else
       names(SE) <- names(mle)     
    
    muhat <- mu(lambda=lamhat, delta=dhat, N=N, eps=eps)
    sig2 <- varTou(lambda=lamhat, delta=dhat, N=N, eps=eps) 
    w <- sig2
	
	   }
    
          
    fit <- list(coefficients = beta,
                 delta = dhat, 
                 se = SE, 
                 vcov = V, 
                 fitted.values = fitval,
                 lambda = lamhat,
                 mu = muhat,
                 var = sig2, 
                 residuals = y - muhat,
                 df = nrow(x) - ncol(x), 
                 logLik= out$value,
                 convergence = out$convergence,
                 start.beta = start.beta, 
                 start.delta = start.delta,
                 w = w,  
                 regress = regress,
                 method=method,
                 data=data)
                 
    class(fit) <- "toureg"
    return(fit)
}


