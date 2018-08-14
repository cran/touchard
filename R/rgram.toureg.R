rgram.toureg <- function (object, xlim = NULL, ylim = NULL, xlab = "Count", ylab = NULL, 
    main = NULL, breaks = NULL, border = "black", fill = "lightgray", col = "blue", 
    lwd = 2, pch = 19, lty = 1,  axes = TRUE, width = NULL, plot = TRUE, 
    style = c("hanging", "standing", "suspended"), scale = c("sqrt", "raw"), max = NULL, ...) 
{
    mt <- terms(object)
    mf <- model.frame(object)
    y <- model.response(mf)
    lhat <- object$fitted
    dhat <- object$delta
    if (is.null(max)){
        max0 <- max(1.5 * max(y), 20L)
    }else{  
        max0 <- max
	}
	
   obs <- as.vector(xtabs( ~ factor(y, levels = 0L:max0)))
   x <- 0L:max0
   p <- matrix(NA, length(lhat), length(x))

   for (i in x){
	 p[, i + 1L] <- dtouch(i, lambda = lhat, delta = dhat)
   }

   exp <- colSums(p)

   scale <- match.arg(scale)
   style <- match.arg(style)
   
    if(is.null(main)){
		main <- if(scale == "sqrt") 
		    paste("Rootogram", " (",style,")", sep="")
        else
		    paste("Histogram", " (",style,")", sep="")
	}
   
   if (is.null(xlab)) 
        xlab <- as.character(attr(mt, "variables"))[2L]
   if (is.null(main)) 
        main <- deparse(substitute(object))
   
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

