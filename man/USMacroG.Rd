\name{USMacroG}
\alias{USMacroG}
\title{US Macroeconomic Data (1950--2000, Greene)}

\description{
Time-series data on 12 US macroeconomic variables for 1950--2000.
}

\usage{data("USMacroG")}

\format{
A quarterly multiple time series from 1950(1) to 2000(4) with 12 variables.
 \describe{
    \item{gdp}{Real gross domestic product (in billion USD),}
    \item{consumption}{Real consumption expenditures,}
    \item{invest}{Real investment by private sector,}
    \item{government}{Real government expenditures,}
    \item{dpi}{Real disposable personal income,}
    \item{cpi}{Consumer price index,}
    \item{m1}{Nominal money stock,}
    \item{tbill}{Quarterly average of month end 90 day treasury bill rate,}
    \item{unemp}{Unemployment rate,}
    \item{population}{Population (in million), interpolation of year end figures
      using constant growth rate per quarter,}
    \item{inflation}{Inflation rate,}
    \item{interest}{Ex post real interest rate (essentially, \code{tbill - inflation}).} 
 }
}

\source{
Online complements to Greene (2003). Table F5.1.

\url{http://pages.stern.nyu.edu/~wgreene/Text/tables/tablelist5.htm}
}

\references{
Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Upper Saddle River, NJ: Prentice Hall.
}

\seealso{\code{\link{Greene2003}}, \code{\link{USMacroSW}}, \code{\link{USMacroSWQ}},
  \code{\link{USMacroSWM}}, \code{\link{USMacroB}}}

\examples{
## data and trend as used by Greene (2003)
data("USMacroG")
USMacroG <- as.ts(merge(as.zoo(USMacroG), trend = 1:nrow(USMacroG) - 1))

## Example 6.1
## Table 6.1
library("dynlm")
fm6.1 <- dynlm(log(invest) ~ tbill + inflation + log(gdp) + trend, data = USMacroG)
fm6.3 <- dynlm(log(invest) ~ I(tbill - inflation) + log(gdp) + trend, data = USMacroG)
summary(fm6.1)
summary(fm6.3)
deviance(fm6.1)
deviance(fm6.3)
vcov(fm6.1)[2,3] 

## F test
linear.hypothesis(fm6.1, "tbill + inflation = 0")
## alternatively
anova(fm6.1, fm6.3)
## t statistic
sqrt(anova(fm6.1, fm6.3)[2,5])
 
## Example 8.2
## Ct = b0 + b1*Yt + b2*Y(t-1) + v
fm1 <- dynlm(consumption ~ dpi + L(dpi), data = USMacroG)
## Ct = a0 + a1*Yt + a2*C(t-1) + u
fm2 <- dynlm(consumption ~ dpi + L(consumption), data = USMacroG)

## Cox test in both directions:
coxtest(fm1, fm2)
## ...and do the same for jtest() and encomptest().
## Notice that in this particular case two of them are coincident.
jtest(fm1, fm2)
encomptest(fm1, fm2)
## encomptest could also be performed `by hand' via
fmE <- dynlm(consumption ~ dpi + L(dpi) + L(consumption), data = USMacroG)
waldtest(fm1, fmE, fm2)

## More examples can be found in:
## help("Greene2003")
}

\keyword{datasets}
