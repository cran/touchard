ptouch0 <- function (x, lambda, delta, N=NULL, eps=sqrt(.Machine$double.eps)) sum(dtouch(0:x, lambda, delta, N, eps))
ptouch <- Vectorize(ptouch0, 'x')

