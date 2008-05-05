\name{DutchAdvert}
\alias{DutchAdvert}
\title{TV and Radio Advertising Expenditures Data}

\description{
Time series of television and radio advertising expenditures (in real terms)
in The Netherlands.
}

\usage{data("DutchAdvert")}

\format{
A four-weekly multiple time series from 1978(1) to 1994(13) with 2 variables.
\describe{
  \item{tv}{Television advertising expenditures.}
  \item{radio}{Radio advertising expenditures.}
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
data("DutchAdvert")
plot(DutchAdvert)

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

## Franses (1998), Table 5.4
round(eacf(log(DutchAdvert[,"tv"]), lag = c(1:19, 26, 39)), digits = 3)
}

\keyword{datasets}
