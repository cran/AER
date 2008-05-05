\name{ArgentinaCPI}
\alias{ArgentinaCPI}
\title{Consumer Price Index in Argentina}
\description{
Time series of consumer price index (CPI) in Argentina (index with 1969(4) = 1).
}
\usage{data("ArgentinaCPI")}

\format{
A quarterly univariate time series from 1970(1) to 1989(4).
}

\source{
Online complements to Franses (1998).

\url{http://www.few.eur.nl/few/people/franses/research/book2.htm}
}

\references{
De Ruyter van Steveninck, M.A. (1996). \emph{The Impact of Capital Imports; Argentina 1970--1989}.
Amsterdam: Thesis Publishers.

Franses, P.H. (1998). \emph{Time Series Models for Business and Economic Forecasting}.
Cambridge, UK: Cambridge University Press.
}

\seealso{\code{\link{Franses1998}}}

\examples{
data("ArgentinaCPI")
plot(ArgentinaCPI)
plot(log(ArgentinaCPI))
}

\keyword{datasets}
