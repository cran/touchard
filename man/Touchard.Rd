\name{Touchard}
\alias{Touchard}
\alias{dtouch}
\alias{ptouch}
\alias{qtouch}
\alias{rtouch}
\alias{tau}

\title{The Touchard Distribution}

\description{Density, normalizing constant, distribution function, quantile function and random number generation for the Touchard distribution 
	with Poisson-like parameter equal to \code{lambda} and shape/dispersion parameter equal to \code{delta}.}

\usage{
dtouch(x, lambda, delta, N=NULL, eps=sqrt(.Machine$double.eps), log = FALSE)
ptouch(x, lambda, delta, N=NULL, eps=sqrt(.Machine$double.eps))
qtouch(p, lambda, delta, N=NULL, eps=sqrt(.Machine$double.eps))
rtouch(n, lambda, delta, N=NULL, eps=sqrt(.Machine$double.eps))
tau(lambda, delta, N=NULL, eps=sqrt(.Machine$double.eps))
}
%- maybe also 'usage' for other objects documented here.
\arguments{

\item{x}{
	vector of quantiles  
}
\item{p}{
	vector of probabilities.}
\item{n}{ 	
number of observations. 
}
\item{lambda}{
	Poisson-like (location) parameter which corresponds to the mean of the distribution when 'delta' = 0  	
}
\item{delta}{
	shape/dispersion parameter which produces unequal dispersion (var x mean) when different from zero and
	mild zero excess compared to the Poisson distribution  
}
\item{N}{ 
	number of terms in the computation (series) of the normalizing constant. If \code{NULL} a recursion formula
	is used and iterated until the specified  relative error is reached.
}

\item{eps}{ relative error in the computation (series) of the normalizing constant. Only used if 'N=NULL'.
	See reference for details.  }

\item{log}{
	logical; if TRUE, probability p is given as log(p).}

}
 

\details{ 
	The Touchard distribution with parameters \eqn{\lambda} (l) and \eqn{\delta} (d) has density
	\deqn{f(x) = \frac{\lambda^x (x+1)^{\delta}}{x! \tau(\lambda,\delta)}}{[ l^x (x+1)^d ] / [ x! tau(l,d) ]}
for \eqn{y=0,1,2,...}, \eqn{\lambda > 0 }  and \eqn{\delta} real. 	


}
\value{
\code{dtouch} gives the density, \code{ptouch} gives the distribution function, \code{qnorm} gives the quantile function, 
and \code{rtouch} generates random deviates.

\code{rtouch} uses the inverse transform method. The length of the result is determined by \code{n} 
and is the maximum of the lengths of the numerical arguments for the other functions.
The numerical arguments other than \code{n} are recycled to the length of the result. 

\code{qtouch} uses an initial approximation based on the Cornish-Fisher expansion followed by
a search in the appropriate direction.
 
\code{tau} gives the value of the normalizing constant in the Touchard density.  

}
\references{
Matsushita RY, Pianto D, Andrade BB, Cancado A, Silva S (2018) The Touchard distribution, 
\emph{Communications in Statistics - Theory and Methods}, <doi:10.1080/03610926.2018.1444177>
}
\author{
\packageAuthor{touchard}
}
%\note{
%%  ~~further notes~~
%}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{rgram}} , \code{\link{touplot}}  
}

\examples{
for(N in c(2, 5, 10, 20, 50)) print( tau(lambda=7, delta=-1, N) )
tau(lambda=7, delta=-1)
dtouch(0:10, lambda=7, delta=-1)
ptouch(0:10, lambda=7, delta=-1)
qtouch(c(.1,.25,.5,.75,.9), lambda=7, delta=-1)
rtouch(10, lambda=7, delta=-1)
}
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
