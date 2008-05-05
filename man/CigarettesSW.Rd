\name{CigarettesSW}
\alias{CigarettesSW}
\title{Cigarette Consumption Panel Data}

\description{
Panel data on cigarette consumption for the 48 continental US States from 1985--1995.
}

\usage{data("CigarettesSW")}

\format{
A data frame containing 48 observations on 7 variables for 2 periods. 
\describe{
\item{state}{Factor indicating state.}
\item{year}{Factor indicating year.}
\item{cpi}{Consumer price index.}
\item{population}{State population.}
\item{packs}{Number of packs per capita.}
\item{income}{State personal income (total, nominal).}
\item{tax}{Average state, federal and average local excise taxes for fiscal year.}
\item{price}{Average price during fiscal year, including sales tax.}
\item{taxs}{Average excise taxes for fiscal year, including sales tax.}
}
}

\source{
Online complements to Stock and Watson (2007). 

\url{http://wps.aw.com/aw_stock_ie_2/}
}

\references{
Stock, J.H. and Watson, M.W. (2007). \emph{Introduction to Econometrics},
2nd ed. Boston: Addison Wesley.
}

\seealso{\code{\link{StockWatson2007}}, \code{\link{CigarettesB}}}

\examples{
## Stock and Watson (2007)
## data and transformations 
data("CigarettesSW")
CigarettesSW$rprice <- with(CigarettesSW, price/cpi)
CigarettesSW$rincome <- with(CigarettesSW, income/population/cpi)
CigarettesSW$tdiff <- with(CigarettesSW, (taxs - tax)/cpi)
c1985 <- subset(CigarettesSW, year == "1985")
c1995 <- subset(CigarettesSW, year == "1995")

## Equations 12.9--12.11
fm_s1 <- lm(log(rprice) ~ tdiff, data = c1995)
coeftest(fm_s1, vcov = sandwich)
fm_s2 <- lm(log(packs) ~ fitted(fm_s1), data = c1995)
fm_ivreg <- ivreg(log(packs) ~ log(rprice) | tdiff, data = c1995)
coeftest(fm_ivreg, vcov = sandwich) ## s.e.s do not match

## Equation 12.15
fm_ivreg2 <- ivreg(log(packs) ~ log(rprice) + log(rincome) | log(rincome) + tdiff, data = c1995)
coeftest(fm_ivreg2, vcov = sandwich) ## s.e.s do not match
## Equation 12.16
fm_ivreg3 <- ivreg(log(packs) ~ log(rprice) + log(rincome) | log(rincome) + tdiff + I(tax/cpi),
  data = c1995)
coeftest(fm_ivreg3, vcov = sandwich)  ## s.e.s do not match

## More examples can be found in:
## help("StockWatson2007")
}

\keyword{datasets}
