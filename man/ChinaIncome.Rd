\name{ChinaIncome}
\alias{ChinaIncome}
\title{Chinese Real National Income Data}
\description{
Time series of real national income in China per section (index with 1952 = 100).
}
\usage{data("ChinaIncome")}

\format{
An annual multiple time series from 1952 to 1988 with 5 variables.
\describe{
  \item{agriculture}{Real national income in agriculture sector.}
  \item{industry}{Real national income in industry sector.}
  \item{construction}{Real national income in construction sector.}
  \item{transport}{Real national income in transport sector.}
  \item{commerce}{Real national income in commerce sector.}
}
}

\source{
Online complements to Franses (1998).

\url{http://www.few.eur.nl/few/people/franses/research/book2.htm}
}

\seealso{\code{\link{Franses1998}}}

\references{
Chow, G.C. (1993). Capital Formation and Economic Growth in China.
\emph{Quarterly Journal of Economics}, \bold{103}, 809--842.

Franses, P.H. (1998). \emph{Time Series Models for Business and Economic Forecasting}.
Cambridge, UK: Cambridge University Press.
}

\examples{
data("ChinaIncome")
plot(ChinaIncome)
}

\keyword{datasets}
