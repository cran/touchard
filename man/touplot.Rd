\title{Touchardness Plot}
\name{touplot}
\alias{touplot}

\keyword{touchard}

\description{
  Touchardness Plot: diagnostic distribution plot for the Touchard model.
}
\usage{
touplot(x, freq = NULL, plot = TRUE, conf.level = 0.95, 
                       ylab = "Count Metameter", xlab = "Count", 
                       main = "Touchardness Plot", ...)
}\arguments{
  \item{x}{either a vector of counts or a 1-way table of frequencies of counts.}
\item{freq}{
	unnecessary (and ignored) if \code{x} is a table or the raw data; must be provided if \code{x} is 
	the unique values of the observed counts, in the same order.
}
\item{plot}{logical. Should the plot be displayed?}
\item{conf.level}{numeric in (0,1). Confidence level for confidence intervals.}
\item{main}{character. Title for the plot.}
\item{xlab}{character. Label for the x axis.}
\item{ylab}{character. Label for the y axis.}
  \item{\dots}{other parameters to be passed through to plotting functions.}
}

\details{
  Plots the number of occurrences (counts) against the
  count metameter of the Touchard distribution.  
  Circles are the observed count metameters and 
  the filled points show the confidence interval (dashed lines) centers.
  Estimate of \eqn{\lambda} based on the fitted line
  is shown on top margin along with the MLE (\code{\link{toufit}}).  
  If the Touchard model fits the data well, the plot should show a straight line
  and the two estimates of \eqn{\lambda} should be close.
  
  
}
\value{
Returns invisibly a data frame containing the counts (\code{y}), 
frequencies (\code{freq}), count metameter (\code{metameter}), the 
CI center (\code{CIcenter}) and the CI margin (\code{CImargin}).
}
\author{
  \packageAuthor{touchard}
}

\references{
  D. C. Hoaglin, F. Mosteller & J. W. Tukey (eds.),
  \emph{Exploring Data Tables, Trends and Shapes}, chapter 9.
  John Wiley & Sons, New York.

  M. Friendly & D. Meyer (2015),
  \emph{Discrete Data Analysis with R}. Chapman and Hall.
}


\seealso{
\code{\link{rgram}}   
}


\examples{
data(Epilepsy)
touplot(seizures)
touplot(seizures, plot=FALSE)
}
