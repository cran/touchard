\name{rgram}
\alias{rgram}

\alias{rgram.toufit}
\alias{rgram.toureg}

\title{Tukey's (hanging, standing, suspended) Rootogram for Assessing a Touchard Model Fit}

\description{
  The Hanging Rootogram is variation of the histogram with the vertical axis showing the square
root of the frequencies hanging from (Touchard) fitted values so that the discrepancies are
visualized against a straight line (the axis) rather than against a curve.
}

\usage{
rgram(object, xlim = NULL, ylim = NULL, xlab = "Count", ylab = NULL, main = NULL, 
breaks = NULL, border = "black", fill = "lightgray", col = "blue", 
         lwd = 2, pch = 19, lty = 1,  axes = TRUE, width = NULL, plot = TRUE, 
         style = c("hanging", "standing", "suspended"), scale = c("sqrt", "raw"), 
         max = NULL, \dots)

}

\arguments{
  \item{object}{either a \code{toufit} or a \code{toureg} object. }
   \item{breaks}{numeric. Breaks for the bar intervals.}
  \item{style}{character. Style of rootogram (see below).}
  \item{scale}{character. Scale of vertical axis, raw frequencies or their
    square roots; \code{style="standing"} with \code{scale="raw"} 
    results in the usual histogram with the fitted curve superimposed.}
  \item{plot}{logical.  Should the plot be displayed?}
  \item{width}{numeric. Widths of the histogram bars.}
 \item{main}{character. Title for the plot.}
\item{xlab}{character. Label for the x axis.}
\item{ylab}{character. Label for the y axis.}
 \item{xlim, ylim, border, fill, col, lwd, pch, lty, axes}{graphical
    parameters. }
 \item{max}{the largest count value used when using the \code{rgram} 
	 with \code{toureg} objects. If \code{NULL}, it sets to
	 \code{max(1.5 * max(count), 20L)}.}
  \item{\dots}{further graphical parameters passed to the plotting
    function.}
}

\details{
The rootogram is a visual tool for comparing the empirical
distribution and fitted values (here from a Touchard model). 
The square-root scale de-emphasizes outlying values and right skewness 
(common for count data). The hanging (from the fitted values) style allows 
discrepancies to be visualized against a straight line (the axis).

The code has been largely based on more general rootogram functions: 
\code{\link[vcd]{rootogram}} in package \pkg{vcd} and 
\code{countreg::rootogram} available through R-forge.
}


\value{
  Returns invisibly a data frame with quantities used in plotting.
}
\author{
	\packageAuthor{touchard}
}

\references{
Friendly M, Meyer D (2015).
  \emph{Discrete Data Analysis with R}. Chapman and Hall.
  
Kleiber C, Zeileis A (2016).
  Visualizing Count Data Regressions Using Rootograms. 
  \emph{The American Statistician}, \bold{70}(3), 296--303.
  \doi{10.1080/00031305.2016.1173590}.

Tukey JW (1977). \emph{Exploratory Data Analysis}. Addison-Wesley, Reading.

}

\seealso{    \code{\link[vcd]{rootogram}},  \code{\link{touplot}}    } 

\examples{
data(Epilepsy)
seiz <- seizures
fm <- toufit( seiz )
rgram(fm)

data(Accidents)
acc <- within(subset(Accidents, subset=Year==2013), {
FriSat <- ifelse(WeekDay > 5, 1, 0)
Spring <- ifelse(Season == "Spring", 1, 0)
})
fmTraff <- toureg( Y ~ FriSat + Spring, data = acc )
rgram(fmTraff)

}

\keyword{rootogram}
