\name{toureg}
\alias{toureg}
%\alias{toureg.fit}
\alias{gleverage}
\alias{cooks.dist}
\alias{hvalues}
\alias{predict}
\alias{predict.toureg}
\alias{plot.toureg}
\alias{residuals.toureg}
\alias{residuals}

\title{Touchard Regression}

\description{Touchard Regression via either maximum likelihood or quasi-likelihood.}


\usage{
toureg(formula, data, x = FALSE, y = FALSE, start.beta, start.delta, 
        parscale = rep.int(1, length(start.beta)+1),  maxit, abstol = -Inf, 
        reltol = 1e-6, etol = 1e-6, gtol = 1e-4,
        N=100, eps=1e-6, dm = 10, regress = c("mu", "lambda"), 
        method = c("BFGS", "CG", "Nelder-Mead", "glm", "qp1", "qp2"),  \dots )

 }
 
 \arguments{
  \item{formula}{an object of class "formula" (or one that can be coerced to that class): 
	  a symbolic description of the model to be fitted.}
  \item{data}{an optional data frame containing the variables in the model. 
	  If not found in data, the variables are taken from environment(formula), 
	  typically the environment from which \code{toureg} is called. }
  \item{start.beta, start.delta}{starting values for the parameters in the linear predictor. If missing, 
	  the default values explained in the reference below are used.}
  \item{parscale, maxit, abstol, reltol}{arguments passed to \code{control} in \code{\link{optim}}. \code{reltol} is alos used as 
	  a relative tolerance for termination of iteratively weighted least squares (\code{method="glm"}) if the relative absolute change 
	  in the function (\strong{f} = loglikelihood) falls below 
	  \code{reltol}*(\code{reltol} + |\strong{f}|).}
  \item{gtol}{iteratively weighted least squares (\code{method="glm"}) stops if largest component (in magnitude) of gradient is less than 
	  this value.}
  \item{etol}{iteratively weighted least squares (\code{method="glm"}) stops when the
          last relative step length is sufficiently small, i.e. below 
          \code{etol}*(\code{etol} + ||\strong{b}||_2), where \strong{b} is current
          state of the minimizer.}
  \item{dm}{non-zero scalar:  with \code{method} set to \code{"glm"} or \code{"qp2"} the estimated \strong{delta} 
	  is obtained by solving some nonlinear equation with root-finding started in the interval 
	  \code{start.delta} +/- \code{dm}. }
  \item{regress}{whether regression is based on log(\code{"mu"} or log(\code{"lambda"}.}
  \item{method}{optimization method for maximization of loglikelihood: 
	  (i) Broyden-Fletcher-Goldfarb-Shanno, Conjugate Gradient or Nelder-Mead as implemented 
	  in \code{\link{optim}};
	  
	  (ii) iteratively weighted least squares (given \strong{delta}) combined with optimization 
	   over \strong{delta} (given the regression coefficients) as in GLM-type models;
	   
	   or 
	   
	  (iii) quasi-Poisson-Touchard (QPT) method with two variants: 
	  \code{"qp1"} assumes variance = \strong{mu-delta} and \code{"qp2"} assumes
	  the exact Touchard variance (see also \code{dm}).}
  \item{N, eps}{arguments passed to \code{\link{tau}}.}
  \item{x, y}{logical values indicating whether the response
          vector and model matrix used in the fitting process should be
          returned as components of the returned value.}
  
  \item{\dots}{not used.}
}
\details{
  Touchard regression with either \eqn{log(\mu)} or \eqn{log(\lambda)} modeled linearly on the predictors as described 
  in Andrade et al (submitted). 
  Estimation can be performed by maximum likelihood via \code{\link[stats]{optim}} with three available methods 
  (\code{'BFGS', 'CG', 'Nelder-Mead'}) and analytical gradients. Default starting values for the coefficients 
  are obtained from Poisson GLM. Default starting value for \eqn{\delta} is obtained by regressing the metameter on the
  sufficient statistics Y and log(Y+1).  Standard errors are obtained from the
  diagonal of inverse of observed Fisher information as reported at the final iteration.
  
  Estimation may also be performed by combination of iteratively weighted least squares and maximization over \eqn{\delta}
  given current estimate of \eqn{\beta}. Details are given Andrade et al (submitted).
  
  Finally, estimation can be performed by Poisson Quasi-MLE (or Poisson pseudo-MLE): the estimator for is \eqn{\beta} 
  is the same as in the Poisson model (which can be thought of as simply a motivation to 
  the first-order condition defining the estimator); the variance is specified independently without 
  restriction of equidispersion. Two specifications are available: (i) a linear specification variance = \strong{mu-delta}
  which corresponds to an approximation to the Touchard variance and (ii) the exact Touchard variance, both allowing
  for under- and over-dispersion.  Details are given Andrade et al (submitted).
  
  Extractor functions for fitted model objects (of class \code{"toureg"}):
  \code{\link[base]{print}}, \code{\link{summary}}, \code{\link{plot}},
   \code{\link[stats]{residuals}}, \code{\link{predict}}, 
  \code{\link{cooks.distance}}, \code{\link{hatvalues}} and \code{\link{gleverage}}.
 
  \code{toureg} returns an object of class \code{"toureg"}, a list with components as described below.
  %\code{toureg.fit} returns an unclassed list with the same components as \code{toureg} except for \code{call}, \code{terms} and \code{formula}.
  }
  
  \value{
  \item{call}{the original function call.}
  \item{coefficients}{named vector of estimated regression coefficients.}
  \item{convergence}{integer code from \code{otim} indicating either successful completion or faulty termination.}
  \item{data}{the data provided in the function call.} 
  \item{delta}{named vector (of length one) of estimated delta parameter.}
  \item{df}{residual degrees of freedom in the fitted model.}
  \item{fitted.values}{a vector of fitted values of lambda.}
  \item{formula}{the formula provided in the function call.}
  \item{lambda}{vector of fitted values of lambda.}
  \item{loglik}{log-likelihood of the fitted model or pseudo-log-likelihood in case os \code{method} 
	  set to \code{"glm"}, \code{"qp1"} or \code{"qp2"}.}
  \item{method}{method used.} 
  \item{mu}{vector of fitted means.}
  \item{residuals}{vector of raw residuals (y - mu).} 
  \item{se}{standard errors of estimated parameters.}
  \item{start.beta, start.delta}{the starting values for the parameters passed to the optimizations routines.}
  \item{w}{weights in the (projection) hat matrix analogous to GLMs.}
  %\item{offset}{a list of offset vectors used (if any),}
  \item{terms}{the 'terms' object used.}
  \item{var}{vector of fitted variances.}
  \item{vcov}{covariance matrix of estimates.}
  \item{x}{if requested, the model matrix.}
  \item{y}{if requested, the response vector.}
  }
  
  \references{
Andrade, BB; Matsushita, RY; Oliveira, SB (submitted) Analyzing Count Data with the Touchard Model in R. 
\emph{available upon request}.
}

\author{
\packageAuthor{touchard}
}

\seealso{
	\code{\link{glm}}, \code{\link{formula}}
}
  
  \examples{

### Horseshoe crab data used by several textbook sources	  
data(Crabs)

### Model Fitting (with different methods) and Plotting
summary( fm <- toureg(y ~ weight + color, data=Crabs) )
# same as
# summary( fm <- toureg(y ~ weight + color, data=Crabs, regress='lambda', method='BFGS') )

# other methods based on log(mu):
# summary( fm2 <- toureg(y ~ weight + color, data=Crabs, regress='mu', method='glm') )
# summary( fm3 <- toureg(y ~ weight + color, data=Crabs, regress='mu', method='qp1') )


plot(fm)
plot(fm , which = 1)
rgram(fm)

### Diagnostics
plot(hvalues(fm))
plot(gleverage(fm))
plot(cooks.dist(fm))


sum(residuals(fm,'response')^2)
sum(residuals(fm,'pearson')^2)
sum(residuals(fm,'deviance')^2)



### Predicted values for 'newdata' ###
 
# Predicted mean values (on the scale of the response variable, i.e. \hat{\mu}):
predict(fm, newdata=data.frame(weight=c(5,6), color=c(2,4)), type="response", se.fit=TRUE)
# Predicted values of lambda:
predict(fm, newdata=data.frame(weight=c(5,6), color=c(2,4)), type="lambda", se.fit=TRUE)
# Predicted values of the linear predictor  x'beta, SEs not yet available:
predict(fm, newdata=data.frame(weight=c(5,6), color=c(2,4)), type="linpred")
# Predicted variances, i.e. \hat{\sigma}^2, SEs not yet available:
predict(fm, newdata=data.frame(weight=c(5,6), color=c(2,4)), type="variance")


}
