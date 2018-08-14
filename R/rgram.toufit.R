rgram.toufit <- function (object, xlim = NULL, ylim = NULL, xlab = "Count", ylab = NULL, 
    main = NULL, breaks = NULL, border = "black", fill = "lightgray", col = "blue", 
    lwd = 2, pch = 19, lty = 1,  axes = TRUE, width = NULL, plot = TRUE, 
    style = c("hanging", "standing", "suspended"), scale = c("sqrt", "raw"), ...) 
{
    if(object$trunc.at.zero || object$rc) 
    stop("not yet implemented for truncated or right censored data")
   
    scale <- match.arg(scale)
    style <- match.arg(style)
   
    if(is.null(main)){
		main <- if(scale == "sqrt") 
		    paste("Rootogram", " (",style,")", sep="")
        else
		    paste("Histogram", " (",style,")", sep="")
	}
   
	 
    mle <- object$fit$est
    
    # zero padding in _obs_ in both cases
    if(is.null(object$data$freq)){    
		Y <- object$data$x
		MAX <- max(Y)
		x <- 0:MAX
		tab <- table(factor(Y, levels=x))
		obs <- as.vector(tab)
		exp <- length(Y) * dtouch(x, mle[1], mle[2])
 	}else{
		X <- object$data$x
		MAX <- max(X)
        obs <- rep(0, MAX+1)
        obs[X+1] <- object$data$freq
        x <- 0:MAX
        exp <- sum(object$data$freq) * dtouch(x, mle[1], mle[2])
	 }
       
    if (is.null(ylab)) {
        ylab <- if (scale == "raw") 
            "Frequency"
        else "sqrt(Frequency)"
    }
    
    breaks <- (head(x, -1L) + tail(x, -1L))/2
    breaks <- c(2 * head(x, 1L) - head(breaks, 1L), breaks, 
            2 * tail(x, 1L) - tail(breaks, 1L))
            
    if (is.null(width))   width <- 0.9      
    
    if (scale == "sqrt") {
        obs <- sqrt(obs)
        exp <- sqrt(exp)
    }
    
    y <- if (style == "hanging") exp - obs  else  0
    height <- if (style == "suspended")  exp - obs  else  obs
  
    width <- diff(breaks) * width
    line <- exp
    
    val <- data.frame(observed = obs, expected = exp, 
        x = x, y = y, width = diff(breaks) * width, height = height)
    
    if(plot){
      xleft <- x - width/2
      xright <- x + width/2
      ybottom <- y
      ytop <- y + height
      if (is.null(xlim)) 
            xlim <- range(c(xleft, xright))
      if (is.null(ylim)) 
            ylim <- range(c(ybottom, ytop, line))
      plot(0, 0, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, 
            ylab = ylab, main = main, axes = FALSE, ...)
      if (axes) {
           axis(1)
           axis(2)
      }
      rect(xleft, ybottom, xright, ytop, border = border, col = fill)
      abline(h = 0, col = border)
      lines(x, line, col = col, pch = pch, type = "b", 
              lty = lty, lwd = lwd)
      invisible(val)
      }else{
		  return(val)
	  }
	  
           
}

