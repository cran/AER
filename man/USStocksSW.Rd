\name{USStocksSW}
\alias{USStocksSW}

\title{Monthly US Stock Returns (1931--2002, Stock \& Watson)}

\description{
Monthly data from 1931--2002 for US stock prices, measured by the
broad-based (NYSE and AMEX) value-weighted index of stock prices as
constructed by the Center for Research in Security Prices (CRSP). 
}

\usage{data("USStocksSW")}

\format{
A monthly multiple time series from 1931(1) to 2002(12) with 2 variables.
\describe{
  \item{returns}{monthly excess returns. The monthly return on stocks
    (in percentage terms) minus the return on a safe asset (in this case:
    US treasury bill). The return on the stocks includes the price changes
    plus any dividends you receive during the month.}
  \item{dividend}{100 times log(dividend yield). (Multiplication by 100
    means the changes are interpreted as percentage points). It is calculated 
    as the dividends over the past 12 months, divided by the price in the current month.}
  }
}

\source{
Online complements to Stock and Watson (2007). 

\url{http://wps.aw.com/aw_stock_ie_2/}
}

\references{
Campbell, J.Y., and Yogo, M. (2006). Efficient Tests of Stock Return Predictability
\emph{Journal of Financial Economics}, \bold{81}, 27--60.

Stock, J.H. and Watson, M.W. (2007). \emph{Introduction to Econometrics}, 2nd ed. Boston: Addison Wesley.
}

\seealso{\code{\link{StockWatson2007}}}

\examples{
data("USStocksSW")
plot(USStocksSW)

## Stock and Watson, p. 540, Table 14.3
library("dynlm")
fm1 <- dynlm(returns ~ L(returns), data = USStocksSW, start = c(1960,1))
coeftest(fm1, vcov = sandwich)
fm2 <- dynlm(returns ~ L(returns, 1:2), data = USStocksSW, start = c(1960,1))
waldtest(fm2, vcov = sandwich)
fm3 <- dynlm(returns ~ L(returns, 1:4), data = USStocksSW, start = c(1960,1))
waldtest(fm3, vcov = sandwich)

## Stock and Watson, p. 574, Table 14.7
fm4 <- dynlm(returns ~ L(returns) + L(d(dividend)), data = USStocksSW, start = c(1960, 1))
fm5 <- dynlm(returns ~ L(returns, 1:2) + L(d(dividend), 1:2), data = USStocksSW, start = c(1960,1))
fm6 <- dynlm(returns ~ L(returns) + L(dividend), data = USStocksSW, start = c(1960,1))
}

\keyword{datasets}
