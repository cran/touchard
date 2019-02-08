\name{score.delta}
\alias{score.delta}

\title{Score Test for Poisson vs Touchard}

\description{
	Score test for the extra (with respect to the Poisson) parameter in the Touchard 
	model (Null: \eqn{\delta} = 0) performed without the need to fit
	the Touchard model. 	Analogous likelihood-ratio and Wald tests (which require fitting
	the full model) are available in the output of \code{\link{toureg}}. 
	}

\usage{
score.delta(x, freq = NULL, max = 50, data)
}

\arguments{
  \item{x}{
	  for univariate data: either an object of class \code{table} with the observed counts as classifying factors \emph{or} the observed
counts (raw data) \emph{or} the unique values of the observed counts in which case the observed
frequencies must be given in argument \code{freq}; for regression data: an object of class \code{"formula"} with the symbolic description 
of the model to be fitted.
}
\item{freq}{
	unnecessary (and ignored) if \code{x} is a table or the raw data or a formula; must be provided if \code{x} is 
	the unique values of the observed counts, in the same order.
}
\item{max}{ 
	number of Poisson terms used in calculations of necessary moments
	}
  \item{data}{data frame containing the variables in the model.}
 }

% \details{}


\value{
  A list with elements
  \item{stat}{the value the chi-squared test statistic.}
  \item{pval}{the p-value for the test (from chi^2(df=1)).}
}
\author{
	\packageAuthor{touchard}
}

\references{
Andrade, BB; Matsushita, RY; Oliveira, SB (submitted) Analyzing Count Data with the Touchard Model. 
\emph{available upon request}.
}

%\seealso{    \code{\link[vcd]{rootogram}},  \code{\link{touplot}}    } 

\examples{
data(Epilepsy)
score.delta(seizures)

data(Crabs)
score.delta( y ~ weight + color, data = Crabs ) 

}

%\keyword{}
