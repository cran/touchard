\name{toureg}
\alias{toureg}
%\alias{toureg.fit}
\alias{gleverage}
\alias{cooks.distance}
\alias{hatvalues}
\alias{predict}
\alias{predict.toureg}
\alias{plot.toureg}


\title{Touchard Regression}

\description{Touchard Regression via maximum likelihood
  with log link for the Poisson parameter (lambda) and precision parameter (delta) constant.
}


\usage{
toureg(formula, data, x = FALSE, y = FALSE, start.beta, start.delta, 
        parscale = rep.int(1, length(start.beta)+1),  maxit = 100L, 
        abstol = -Inf, reltol = sqrt(.Machine$double.eps), N=100, eps=1e-8, 
         \dots)

 }
 
 \arguments{
  \item{formula}{an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
  \item{data}{an optional data frame containing the variables in the model. If not found in data, the variables are taken from environment(formula), 
	  typically the environment from which \code{toureg} is called. }
  \item{start.beta, start.delta}{starting values for the parameters in the linear predictor. If missing, the default values explained in the reference 
are used.}
\item{parscale, maxit, abstol, reltol}{arguments passed to \code{control} in \code{\link{optim}}.}
  \item{N, eps}{arguments passed to \code{\link{tau}}.}
  \item{x, y}{logical values indicating whether the response
          vector and model matrix used in the fitting process should be
          returned as components of the returned value.}
  \item{\dots}{not used.}
}
\details{
  Touchard regression with \eqn{\log(\lambda)} modeled linearly on the predictors as described 
  in Andrade et al (submitted). \eqn{\delta} is, so far, kept constant. 
  Estimation is performed by maximum likelihood via \code{\link[stats]{optim}} using
  \code{method = 'BFGS'} and analytical gradients. Default starting values for the coefficients 
  are obtained from Poisson GLM. Default starting value for \eqn{\delta} is obtained by regressing the metameter on the
  sufficient statistics Y and log(Y+1).  Standard errors are obtained from the
  diagonal of inverse of observed Fisher information as reported at the final iteration.
  
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
  \item{loglik}{log-likelihood of the fitted model.}
  \item{mu}{vector of fitted means.}
  \item{residuals}{vector of pearson residuals (y - mu)/sigma.} 
  \item{se}{standard errors of estimated regression coefficients.}
  \item{start.beta, start.delta}{the starting values for the parameters passed to the \code{optim} call.}
  %\item{weights}{the weights used (if any),}
  %\item{offset}{a list of offset vectors used (if any),}
  \item{terms}{the 'terms' object used.}
  \item{var}{vector of fitted varainces.}
  \item{vcov}{covariance matrix of all parameters (beta and delta) in the model.}
  \item{x}{if requested, the model matrix.}
  \item{y}{if requested, the response vector.}
  }
  
  \references{
Andrade, BB; Matsushita, RY; Oliveira, SB (submitted) Analyzing Count Data with the Touchard Model in R. 
\emph{J. of Statistical Software}.
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

### Model Fitting and Plotting
summary( fm <- toureg(y ~ weight + color, data=Crabs) )

plot(fm)
plot(fm , which = 1)
rgram(fm)

### Diagnostics
plot(hatvalues(fm))
plot(gleverage(fm))
plot(cooks.distance(fm))

### Predicted values for 'newdata' ###
 
# The default _type="invlink"_  is on the scale of the inverse link, i.e. \hat{lambda}=exp(x*beta); 
predict(fm, newdata=data.frame(weight=c(5,6), color=c(2,4)), se.fit=TRUE)
# the alternative "response" is on the scale of the response variable, i.e. \hat{mu}, 
# which is a function of \hat{lambda} and \hat{delta};
predict(fm, newdata=data.frame(weight=c(5,6), color=c(2,4)), type="response", se.fit=TRUE)
# and finally the alternative "variance" yields predicted variances, i.e. \hat{sigma}^2,  
# which is a function of \hat{lambda} and \hat{delta} and for which SEs are not available so far.
predict(fm, newdata=data.frame(weight=c(5,6), color=c(2,4)), type="variance")


}
