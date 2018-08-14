touplot <- function(x, freq = NULL, plot = TRUE, conf.level = 0.95, 
                       ylab = "Count Metameter", xlab = "Count", 
                       main = "Touchardness Plot", ...){
   if(conf.level<0 || conf.level>1) stop("_conf.level_ not in (0,1)") 
   if(is.table(x)){
	 if(!is.null(freq)) warning("_freq_ has been ignored since _x_ is already of class 'table'")
	 ny <- as.vector(x)
	 y <- as.numeric(names(x))
     }else{
	   if(is.null(freq)){ 
		   X <- table(x)
		   ny <- as.vector(X)
	       y <- as.numeric(names(X))
	       }else{
			  ny <- freq  
			  y <- x
		   }
		   }

    n <- sum(ny)
    mle <- toufit(x)$fit$estimate
    del <- mle['delta'][[1]]
    phi <- lgamma(y+1) + log(ny) - log(n) - del*log(y+1)
    nystar <- ifelse(ny > 1.5, ny - 0.67 - .8*ny/n, 1/exp(1))
    phistar <- lgamma(y+1) + log(nystar) - log(n) - del*log(y+1)
    phat <- ny/n
    e <-  sqrt(1 - phat)/sqrt(ny - (0.25 * phat + 0.47) * sqrt(ny))
    qz <- qnorm(1 - (1 - conf.level)/2)
    h <- qz*e
    val <- data.frame(y=y, freq=ny, metameter=phi, CIcenter=phistar, CImargin=h) 
    
    myline <- coef(lm(phi~y))
    b <- myline[2]
    
    l1 <- round(exp(b),2)
    l2 <- round(mle['lambda'][[1]],2)
    
    if(plot){
       plotrix::plotCI(x=y, y=phistar, uiw=h, 
              sfrac=0.01, slty=3, pch=19, pt.bg=par("bg"), 
              ylab=ylab, xlab=xlab, main=main)
       mtext(bquote(
             paste(hat(lambda), ": exp(slope) = ", .(l1), "  |  mle = ", .(l2))
            ), side = 3)
       abline(myline, lwd=2, col="blue")
       points(x=y, y=phi, pch=21)
       invisible(val)
    }else{
       return(val)
   }
}
