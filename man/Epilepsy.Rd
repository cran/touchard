\name{Epilepsy}
\alias{Epilepsy}
\alias{seizures}
\docType{data}
\title{Epilepsy Data}
\description{Data on Epileptic Seizures.}
\usage{data("Epilepsy")}
\format{
  A data frame with 351 observations on counts of epileptic seizures.
  \describe{
\item{\code{seizures}}{numeric vector of length 351.}
  }
}
%%\details{}
\source{
see reference below
}
\references{
Bhati D, Sastry DVS, Qadri PZM (2015). A New Generalized Poisson-Lindley Distribution:
Applications and Properties. \emph{Austrian Journal of Statistics}, 4, 35-51.
}
\examples{
data(Epilepsy)
seizures
table(seizures)
}
\keyword{datasets}
