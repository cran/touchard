\name{toufit}
\alias{toufit}
\title{Touchard Estimation}

\description{Maximum-likelihood and method-of-moments estimation of the Touchard model.}

\usage{
toufit(x, freq = NULL, start, method = c("ml", "mm", "gmm"), 
    rc = FALSE, trunc.at.zero = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

\item{x}{
either an object of class \code{table} with the observed counts as classifying factors \bold{or} the observed
counts (raw data) \bold{or} the unique values of the observed counts in which case the observed
frequencies must be given in argument \code{freq}; see example below for the three formats.  
}
\item{freq}{
	unnecessary (and ignored) if \code{x} is a table or the raw data; must be provided if \code{x} is 
	the unique values of the observed counts, in the same order.
}
\item{start}{
starting values to be used by optimization/nonlinear solving algorithms. If missing, the default value is
explained in the reference below.
}
\item{method}{
one of 'ml' (maximum likelihood), 'mm' (method of moments based on the usual first and second moment conditions)
or 'gmm' (generalized method of moments based on the first two moments and an extra condition on the first moment
of log(Y+1). 
}
\item{rc}{
TRUE for right-censored data. Not implemented for methods 'mm' and 'gmm'.
}

\item{trunc.at.zero}{ TRUE for model truncated at zero. Not implemented for methods 'mm' and 'gmm'. }


}
 

\details{ For method 'ml' the \code{\link{Touchard}} likelihood is numerically maximized via \code{optim()}
	 with 'method=L-BFGS-B'. For methods 'mm' and 'gmm', the system of moment conditions
	 is solved by calling  \code{nleqslv} and \code{numDeriv}. 
}
\value{
	A list with the following:
	\item{fit}{a list with point estimates, standard errors and variance matrix of estimates}
	\item{aic}{Akaike's information criterion, -2*log-likelihood + 2*p} 
	\item{bic}{Schwarz's Bayesian criterion, -2*log-likelihood + log(n)*p}
	\item{test}{a data frame with likelihood ratio and Wald test results for the Null: delta = 0 (Poisson)}
	\item{method}{estimation method used ('ml', 'mm' or 'gmm')} 
	\item{data}{a list with elements \code{x} and \code{freq}} 
}
\references{
Matsushita RY, Pianto D, Andrade BB, Cancado A, Silva S (2018) The Touchard distribution, 
\emph{Communications in Statistics - Theory and Methods}, <doi:10.1080/03610926.2018.1444177>
 
Andrade, BB; Matsushita, RY; Oliveira, SB (submitted) Analyzing Count Data with the Touchard Model in R. 
\emph{J. of Statistical Software}.
}
\author{
\packageAuthor{touchard}
}
%\note{
%%  ~~further notes~~
%}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{optim}} , \code{\link{rgram}} , \code{\link{Touchard}} , \code{\link{touplot}}  
}

\examples{
data(Epilepsy)
# 'seizures' is vector of raw data (351 counts) 
toufit( x = seizures, method = 'ml' )
toufit( x = seizures, method = 'mm' )
toufit( x = seizures, method = 'gmm' )

# suppose data were recorded as table object
TAB <- table(seizures)
TAB
toufit(TAB)  # same as above (method = 'ml') 

# suppose data were recorded as data.frame 
DF <- data.frame( y = as.numeric(names(TAB)), fr = as.numeric(TAB) )
DF
toufit(x = DF$y, freq = DF$fr) # same as above (method = 'ml') 
}
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
