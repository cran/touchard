\name{Accidents}
\alias{Accidents}
\docType{data}
\title{Trafic Accidents Data}
\description{Data on daily traffic accidents.}
\usage{data("Accidents")}
\format{
  A data frame with 1096 observations on 6 variables.
  \describe{
\item{\code{Year}}{year (2011, 12, 13).}
\item{\code{Month}}{integer, 1--12.}
\item{\code{WeekDay}}{integer, 1--7.}
\item{\code{Season}}{factor (Fall, Spring, Summer, Winter).}
\item{\code{Holiday}}{integer (0=no, 1=yes).}
\item{\code{Y}}{integer, observed count of accidents.}
  }
}
%%\details{}
\source{
<https://data.ny.gov>   Open Data NY: State of New York 
}
%%\references{}
\examples{
data(Accidents)
summary(Accidents)
}
\keyword{datasets}
