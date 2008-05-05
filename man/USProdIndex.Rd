\name{USProdIndex}
\alias{USProdIndex}
\title{Index of US Industrial Production}
\description{
Index of US industrial production (1985 = 100).
}
\usage{data("USProdIndex")}

\format{
A quarterly multiple time series from 1960(1) to 1981(4) with 2 variables.
\describe{
  \item{unadjusted}{raw index of industrial production,}
  \item{adjusted}{seasonally adjusted index.}
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
data("USProdIndex")
plot(USProdIndex, plot.type = "single", col = 1:2)

## EACF tables (Franses 1998, p. 99)
ctrafo <- function(x) residuals(lm(x ~ factor(cycle(x))))
ddiff <- function(x) diff(diff(x, frequency(x)), 1)
eacf <- function(y, lag = 12) {
  stopifnot(all(lag > 0))
  if(length(lag) < 2) lag <- 1:lag
  rval <- sapply(
    list(y = y, dy = diff(y), cdy = ctrafo(diff(y)),
         Dy = diff(y, frequency(y)), dDy = ddiff(y)),
    function(x) acf(x, plot = FALSE, lag.max = max(lag))$acf[lag + 1])
  rownames(rval) <- lag
  return(rval)
}

## Franses (1998), Table 5.1
round(eacf(log(USProdIndex[,1])), digits = 3)

## Franses (1998), Equation 5.6: Unrestricted airline model
## (Franses: ma1 = 0.388 (0.063), ma4 = -0.739 (0.060), ma5 = -0.452 (0.069))
arima(log(USProdIndex[,1]), c(0, 1, 5), c(0, 1, 0), fixed = c(NA, 0, 0, NA, NA))
}

\keyword{datasets}
