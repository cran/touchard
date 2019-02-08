summary.toureg <- function(object, ...)
{
  k <- length(coef(object))	
  se <- object$se[1:k]
  zval <- coef(object) / se
  TAB <- cbind(Estimate = coef(object),
               StdErr = se,
               z.value = zval,
               p.value = 2*pnorm(-abs(zval),mean=0,sd=1))
  res <- list(call=object$call, coefficients=TAB, se=object$se, delta=object$delta, logLik=object$logLik,
              AIC=-2*object$logLik + 2*k, BIC = -2*object$logLik + log(length(object$fitted))*k,
              method=object$method
              #Deviance=object$deviance
  )
  class(res) <- "summary.toureg"
  res
}
