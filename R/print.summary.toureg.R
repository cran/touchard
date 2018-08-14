print.summary.toureg <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")

  printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
  cat("\n")
  
  cat("delta = ", x$delta, "\n")
  cat("SE(delta) = ", x$se['delta'], "\n")
  
  cat("Log-likelihood = ", x$logLik, "\n")
  cat("AIC = ", x$AIC, "   ",  "BIC = ", x$BIC, "\n")

}
