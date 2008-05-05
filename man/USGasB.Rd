\name{USGasB}
\alias{USGasB}
\title{US Gasoline Market Data (1950--1987, Baltagi)}
\description{
Time-series data on the US gasoline market.
}
\usage{data("USGasB")}
\format{
An annual multiple time series from 1950 to 1987 with 6 variables.
\describe{
  \item{cars}{Stock of cars.}
  \item{gas}{Consumption of motor gasoline (in 1000 gallons).}
  \item{price}{Retail price of motor gasoline.}
  \item{population}{Population.}
  \item{gnp}{Real gross national product (in 1982 dollars).}
  \item{deflator}{GNP deflator (1982 = 100).}
}
}

\source{
The data is from Baltagi (2002) and available at

\url{http://www.springeronline.com/sgw/cda/frontpage/0,10735,4-165-2-107420-0,00.html}
}

\references{
Baltagi, B.H. (2002). \emph{Econometrics}, 3rd ed. Berlin, Springer. 

}

\seealso{\code{\link{Baltagi2002}}, \code{\link{USGasG}}}

\examples{
data("USGasB")
plot(USGasB)
}

\keyword{datasets}
