\name{ConsumerGood}
\alias{ConsumerGood}
\title{Properties of a Fast-Moving Consumer Good}
\description{
Time series of distribution, market share and price of a fast-moving
consumer good.
}
\usage{data("ConsumerGood")}

\format{
A weekly multiple time series from 1989(11) to 1991(9) with 3 variables.
\describe{
  \item{distribution}{Distribution.}
  \item{share}{Market share.}
  \item{price}{Price.}
}
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
data("ConsumerGood")
plot(ConsumerGood)
}

\keyword{datasets}
