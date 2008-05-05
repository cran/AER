\name{GermanUnemployment}
\alias{GermanUnemployment}
\title{Unemployment in Germany Data}

\description{
Time series of unemployment rate (in percent) in Germany.
}

\usage{data("GermanUnemployment")}

\format{
A quarterly multiple time series from 1962(1) to 1991(4) with 2 variables.
\describe{
  \item{unadjusted}{Raw unemployment rate,}
  \item{adjusted}{Seasonally adjusted rate.}
}
}

\source{
Online complements to Franses (1998).

\url{http://www.few.eur.nl/few/people/franses/research/book2.htm}
}

\seealso{\code{\link{Franses1998}}}

\references{
Franses, P.H. (1998). \emph{Time Series Models for Business and Economic Forecasting}.
Cambridge, UK: Cambridge University Press.
}

\examples{
data("GermanUnemployment")
plot(GermanUnemployment, plot.type = "single", col = 1:2)
}

\keyword{datasets}
