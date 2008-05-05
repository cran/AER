\name{DJFranses}
\alias{DJFranses}
\title{Dow Jones Index Data (Franses)}
\description{
Dow Jones index time series computed at the end of the week where
week is assumed to run from Thursday to Wednesday.
}
\usage{data("DJFranses")}

\format{
A weekly univariate time series from 1980(1) to 1994(42).
}

\source{
Online complements to Franses (1998).

\url{http://www.few.eur.nl/few/people/franses/research/book2.htm}
}

\references{
Franses, P.H. (1998). \emph{Time Series Models for Business and Economic Forecasting}.
Cambridge, UK: Cambridge University Press.
}

\seealso{\code{\link{Franses1998}}}

\examples{
data("DJFranses")
plot(DJFranses)
}

\keyword{datasets}
