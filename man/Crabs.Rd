\name{Crabs}
\alias{Crabs}
\docType{data}
\title{Crabs Data}
\description{Horseshoe crab data in Table 4.3.2 (Agresti, Categ. Data Anal., 3nd edition) by courtesy of Jane Brockmann, Zoology Department, University of Florida.}
\usage{data("Crabs")}
\format{
  A data frame with 173 observations on the following 6 variables.
  \describe{
\item{\code{crab}}{crab id}
\item{\code{y}}{number of satellites}
\item{\code{weight}}{in kg}
\item{\code{width}}{carapace width in cm}		
\item{\code{color}}{has values 1-4 with 1=light}
\item{\code{spine}}{spine condition}
  }
}
%%\details{}
\source{
<http://www.stat.ufl.edu/~aa/cda/data.html>
}
\references{
Agresti, A. (2013). \emph{Categorical Data Analysis}. Wiley, 3d ed.

Brockmann, H. J. (1996), Satellite Male Groups in Horseshoe Crabs, Limulus polyphemus. \emph{Ethology}, 102: 1-21. 

}
\examples{
data(Crabs)
dim(Crabs)
head(Crabs)

}
\keyword{datasets}
