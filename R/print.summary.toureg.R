print.summary.toureg <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")

  printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
  cat("\n")
  
  cat("delta = ", x$delta, "\n")
  if(x$method=="qp1" || x$method=="qp2"){
    cat("SE(delta) = NULL", "\n")
    cat("(Pseudo)Log-likelihood = ", x$logLik, "\n")
    cat("(Pseudo)AIC = ", x$AIC, "   ",  "(Pseudo)BIC = ", x$BIC, "\n")
  }else{
    cat("SE(delta) = ", x$se['delta'], "\n")
    cat("Log-likelihood = ", x$logLik, "\n")
    cat("AIC = ", x$AIC, "   ",  "BIC = ", x$BIC, "\n")
   }

}
