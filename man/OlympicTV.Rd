\name{OlympicTV}
\alias{OlympicTV}
\title{Television Rights for Olympic Games}
\description{
Television rights for Olympic Games for US networks (in millions USD).
}
\usage{data("OlympicTV")}

\format{
A data frame with 10 observations and 2 variables.
\describe{
  \item{rights}{time series of television rights (in million USD),}
  \item{network}{factor coding television network.}
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
data("OlympicTV")
plot(OlympicTV$rights)
}

\keyword{datasets}
