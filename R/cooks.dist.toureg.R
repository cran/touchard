
cooks.dist.toureg <- function(object, ...){
	h <- gleverage(object)   # former hatvalues
	rp <- residuals(object)
	k <- length(object$coefficients)
	cd <- (h*rp^2) / ( k * (1-h)^2 )
	return(cd)
}
