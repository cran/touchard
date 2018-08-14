print.toureg <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n Coefficients:\n")
  print(x$coefficients)
  cat("\n delta:", x$delta, "\n")
}
