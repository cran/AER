\name{USConsump1979}
\alias{USConsump1979}
\title{US Consumption Data (1970--1979)}
\description{
Time-series data on US income and consumption expenditure, 1970--1979.
}
\usage{data("USConsump1979")}
\format{
An annual multiple time series from 1970 to 1979 with 2 variables.
\describe{
  \item{income}{Disposable income.} 
  \item{expenditure}{Consumption expenditure.}
}
}

\source{
Online complements to Greene (2003). Table F1.1.

\url{http://pages.stern.nyu.edu/~wgreene/Text/tables/tablelist5.htm}
}

\references{
Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Upper Saddle River, NJ: Prentice Hall.
}

\seealso{\code{\link{Greene2003}}, \code{\link{USConsump1950}}, \code{\link{USConsump1993}}}

\examples{
data("USConsump1979")
plot(USConsump1979)

## Example 1.1 in Greene (2003)
plot(expenditure ~ income, data = as.data.frame(USConsump1979), pch = 19)
fm <- lm(expenditure ~ income, data = as.data.frame(USConsump1979))
summary(fm)
abline(fm)
}

\keyword{datasets}
