\name{USConsump1993}
\alias{USConsump1993}
\title{US Consumption Data (1950--1993)}
\description{
Time-series data on US income and consumption expenditure, 1950--1993.
}
\usage{data("USConsump1993")}

\format{
An annual multiple time series from 1950 to 1993 with 2 variables.
\describe{
  \item{income}{Disposable personal income (in 1987 USD).}
  \item{expenditure}{Personal consumption expenditures (in 1987 USD).}
}
}

\source{
The data is from Baltagi (2002) and available at

\url{http://www.springeronline.com/sgw/cda/frontpage/0,10735,4-165-2-107420-0,00.html}
}

\references{
Baltagi, B.H. (2002). \emph{Econometrics}, 3rd ed. Berlin, Springer. 
}

\seealso{\code{\link{Baltagi2002}}, \code{\link{USConsump1950}}, \code{\link{USConsump1979}}}

\examples{
## data from Baltagi (2002)
data("USConsump1993")
plot(USConsump1993, plot.type = "single", col = 1:2)

## Chapter 5 (p. 122-125)
fm <- lm(expenditure ~ income, data = USConsump1993)
summary(fm)
## Durbin-Watson test (p. 122)
dwtest(fm)
## Breusch-Godfrey test (Table 5.4, p. 124)
bgtest(fm)
## Newey-West standard errors (Table 5.5, p. 125)
coeftest(fm, vcov = NeweyWest(fm, lag = 3, prewhite = FALSE, adjust = TRUE)) 

## Chapter 8
## Recursive residuals
rr <- recresid(fm)
rr
## Recursive CUSUM test
rcus <- efp(expenditure ~ income, data = USConsump1993)
plot(rcus)
sctest(rcus)
## Harvey-Collier test
harvtest(fm)
## NOTE" Mistake in Baltagi (2002) who computes
## the t-statistic incorrectly as 0.0733 via
mean(rr)/sd(rr)/sqrt(length(rr))
## whereas it should be (as in harvtest)
mean(rr)/sd(rr) * sqrt(length(rr))

## Rainbow test
raintest(fm, center = 23)

## J test for non-nested models
library("dynlm")
fm1 <- dynlm(expenditure ~ income + L(income), data = USConsump1993)
fm2 <- dynlm(expenditure ~ income + L(expenditure), data = USConsump1993)
jtest(fm1, fm2)

## Chapter 14
## ACF and PACF for expenditures and first differences
exps <- USConsump1993[, "expenditure"]
(acf(exps))
(pacf(exps))
(acf(diff(exps)))
(pacf(diff(exps)))

## dynamic regressions, eq. (14.8)
fm <- dynlm(d(exps) ~ I(time(exps) - 1949) + L(exps))
summary(fm)
}

\keyword{datasets}
