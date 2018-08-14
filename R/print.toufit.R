print.toufit <- function(x, ...) {
	fm <- x
 if(fm$method == 'ml'){
	cat('   Maximum Likelihood', '\n')
	cat("$fit", '\n')
	print(fm$fit)
	
	cat("AIC = ", fm$aic, "   ",  "BIC = ", fm$bic, "\n")

	cat('\n', 'Null: delta = 0 (Poisson)', '\n')
	printCoefmat(fm$test, P.values=TRUE, has.Pvalue=TRUE)
	}
	
 if( fm$method %in% c('mm','gmm') ){
	r <- ifelse( fm$method == 'mm', 2, 3 ) 
	cat('   Method of Moments', ' ( r = ', r, ')', '\n')
	cat("$fit", '\n')
	print(fm$fit)
	
	cat('\n', 'Null: delta = 0 (Poisson)', '\n')
	#cat('        (Wald Test)', '\n')
	printCoefmat(fm$test, P.values=TRUE, has.Pvalue=TRUE)
	}
 
}
