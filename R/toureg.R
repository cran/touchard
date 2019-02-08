toureg <- function(formula, data, x = FALSE, y = FALSE, start.beta, start.delta, parscale = rep.int(1, length(start.beta)+1),  maxit, abstol = -Inf, 
        reltol = 1e-6, etol= 1e-6, gtol = 1e-4, N=100, eps=1e-6, dm = 10, regress = c("mu", "lambda"),  
        method = c("BFGS", "CG", "Nelder-Mead", "glm", "qp1", "qp2"),...) 
{ 
			UseMethod("toureg")
}
