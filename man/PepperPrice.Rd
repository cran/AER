\name{PepperPrice}
\alias{PepperPrice}
\title{Black and White Pepper Prices}
\description{
Time series of average monthly European spot prices for black and
white pepper (fair average quality) in US dollars per ton.
}
\usage{data("PepperPrice")}

\format{
A monthly multiple time series from 1973(10) to 1996(4) with 2 variables.
\describe{
  \item{black}{spot price for black pepper,}
  \item{white}{spot price for white pepper.}
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

\examples{
## data
data("PepperPrice")
plot(PepperPrice, plot.type = "single", col = 1:2)

## package
library("tseries")
library("urca")

## unit root tests
adf.test(log(PepperPrice[, "white"]))
adf.test(diff(log(PepperPrice[, "white"])))
pp.test(log(PepperPrice[, "white"]), type = "Z(t_alpha)")
pepper_ers <- ur.ers(log(PepperPrice[, "white"]),
  type = "DF-GLS", model = "const", lag.max = 4)
summary(pepper_ers)

## stationarity tests
kpss.test(log(PepperPrice[, "white"]))

## cointegration
po.test(log(PepperPrice))
pepper_jo <- ca.jo(log(PepperPrice), ecdet = "const", type = "trace")
summary(pepper_jo)
pepper_jo2 <- ca.jo(log(PepperPrice), ecdet = "const", type = "eigen")
summary(pepper_jo2)
}

\keyword{datasets}
