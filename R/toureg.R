toureg <- function(formula, data, x = FALSE, y = FALSE, start.beta, start.delta, parscale = rep.int(1, length(start.beta)+1),  maxit = 100L, abstol = -Inf, 
        reltol = sqrt(.Machine$double.eps), N=100, eps=1e-8,  ...) 
{ 
			UseMethod("toureg")
}
