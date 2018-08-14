toureg.formula <- function(formula, data=list(), x = FALSE, y = FALSE, ...)
{
  data <- na.omit(data)
  mf <- model.frame(formula=formula, data=data)
  terms <- attr(mf, "terms")
    
  fit <- toureg.default(x=model.matrix(terms, mf), y=model.response(mf), ...)  
  
  fit$terms <- terms
  fit$formula <- formula
  fit$data <- data
  fit$call <- match.call()
  if(y) fit$y <- model.response(mf)
  if(x) fit$x <- model.matrix(terms, mf)
  fit
}
