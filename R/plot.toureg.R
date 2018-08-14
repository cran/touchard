plot.toureg <- function (x, which = 1:4,  main = "", 
              caption = c("Residuals vs linear predictor", 
                          "Normal Q-Q plot of residuals",
                          "Cook's distance", 
                          "Leverage vs predicted values"), 
              sub.caption = paste(deparse(x$call), collapse = "\n"), 
              type = "Pearson", nsim = 50, level = 0.95, envelope = FALSE, ...) 
{       
	  if (!is.numeric(which) || any(which < 1) || any(which > 4))  stop("_which_ must be in 1:4")
	  #x <- object
	  res <- residuals(x)   
      
      n <- length(res)
      #k <- length(x$coefficients)
      
      show <- rep(FALSE, 4)
      show[which] <- TRUE
      
      old.par <- par(no.readonly = TRUE) # all par settings which
                                           # could be changed.
      on.exit(par(old.par))
      
      if(length(which)>1){
        par(mar=c(5.1,4.1,4.1,2.1)+.5)
        par(mfrow=c(2,2))
      }
      
  
      if(show[1]){
          plot(predict(x, type = "invlink")$fit, res, xlab = "Linear predictor", 
            ylab = paste(type, " residual"),  ...)
        mtext(caption[1], 3, 0.25, cex=.8)
        abline(h = 0, lty = 3, col = "gray")
      }
      
	  if(show[2]){		  
        qn <- qnormal(x, nsim = nsim, level = level, envelope = envelope)  # , type = type
        YLIM <- range(qn)
        qqnorm(qn[,1], xlab="Theoretical Quantiles",
                       ylab="Sample Quantiles", ylim=YLIM, pch=19, cex=0.5, main="") 
        if(envelope){
			par(new=TRUE)
			qqnorm(qn[,3],axes=F,xlab="",ylab="",type="l", ylim=YLIM,lty=1, main="")
			par(new=TRUE)
			qqnorm(qn[,4],axes=F,xlab="",ylab="", type="l", ylim=YLIM,lty=1, main="")
			par(new=TRUE)
			qqnorm(qn[,2],axes=F,xlab="", ylab="", type="l", ylim=YLIM,lty=2, main="")
		}
        mtext(caption[2], 3, 0.25, cex=.8) 
	  }
      
      if(show[3]){
		plot(1:n, cooks.distance(x), xlab = "Obs. number", 
		    ylab = "Cook's distance", 
            type = "h", ...)
        mtext(caption[3], 3, 0.25, cex=.8)  
		  }
		  
	  if(show[4]){
		plot(fitted(x), gleverage(x), xlab = "Fitted value", 
            ylab = "Generalized leverage", ...)
        mtext(caption[4], 3, 0.25, cex=.8)
      }
      invisible()

}


