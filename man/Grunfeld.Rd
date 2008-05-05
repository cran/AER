\name{Grunfeld}
\alias{Grunfeld}

\title{Grunfeld's Investment Data}

\description{
Panel data on 11 large US manufacturing firms over 20 years, for the years 1935--1954.
}

\usage{data("Grunfeld")}

\format{
A data frame containing 20 annual observations on 3 variables for 11 firms. 
\describe{
  \item{invest}{Annual real gross investment.}
  \item{value}{Real value of the firm (shares outstanding).}
  \item{capital}{Real value of the capital stock.}
  \item{firm}{factor with 11 levels: \code{"General Motors"}, \code{"US Steel"}, \code{"General Electric"}, 
              \code{"Chrysler"}, \code{"Atlantic Refining"}, \code{"IBM"}, \code{"Union Oil"}, \code{"Westinghouse"},
              \code{"Goodyear"}, \code{"Diamond Match"}, \code{"American Steel"}.}
  \item{year}{Year.}
}
}

\details{
This is a popular data set for teaching purposes. Unfortunately, there exist several 
different versions (see \url{http://www.stanford.edu/~clint/bench/grunfeld.htm} for details). 
In particular, the variants given by Baltagi (2002) and Greene (2003) differ: 
for \code{"US Steel"} (firm 2), Greene has investment for the year 1940 as 261.6 while Baltagi has 361.6, 
and capital for 1946 as 132.6 while Baltagi has 232.6. 

Here, we provide the original data from Grunfeld (1958, Table 9). The data for
the first 10 firms is identical to that of Baltagi (2002), now also used by
Greene (2008).
}

\source{
The data are taken from Grunfeld (1958, Table 9).

The data for the first 10 firms is identical to that of Baltagi (2002) available at

\url{http://www.springeronline.com/sgw/cda/frontpage/0,10735,4-165-2-107420-0,00.html}
}

\references{
Baltagi, B.H. (2002). \emph{Econometrics}, 3rd ed. Berlin, Springer. 

Baltagi, B.H. (2005). \emph{Econometric Analysis of Panel Data}, 3rd ed. Chichester, UK: John Wiley.

Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Upper Saddle River, NJ: Prentice Hall.

Greene, W.H. (2008). \emph{Econometric Analysis}, 6th edition. Upper Saddle River, NJ: Prentice Hall.

Grunfeld, Y. (1958). \emph{The Determinants of Corporate Investment}. 
Unpublished Ph.D. Dissertation, University of Chicago.
}

\seealso{\code{\link{Baltagi2002}}, \code{\link{Greene2003}}}

\examples{
data("Grunfeld")

## Greene (2003)
## subset of data with mistakes
ggr <- subset(Grunfeld, firm \%in\% c("General Motors", "US Steel",
  "General Electric", "Chrysler", "Westinghouse"))
ggr[c(26, 38), 1] <- c(261.6, 645.2)
ggr[32, 3] <- 232.6

## Tab. 14.2, col. "GM"
fm_gm <- lm(invest ~ value + capital, data = ggr, subset = firm == "General Motors")
mean(residuals(fm_gm)^2)   ## Greene uses MLE

## Tab. 14.2, col. "Pooled"
fm_pool <- lm(invest ~ value + capital, data = ggr)

## equivalently
library("plm")
pggr <- plm.data(ggr, c("firm", "year"))
library("systemfit")
fm_ols <- systemfit(invest ~ value + capital, data = pggr, method = "OLS")
fm_pols <- systemfit(invest ~ value + capital, data = pggr, method = "OLS", pooled = TRUE)

## Tab. 14.1
fm_sur <- systemfit(invest ~ value + capital, data = pggr, method = "SUR", methodResidCov = "noDfCor")
fm_psur <- systemfit(invest ~ value + capital, data = pggr, method = "SUR", pooled = TRUE, 
  methodResidCov = "noDfCor", residCovWeighted = TRUE)

## More examples can be found in:
## help("Greene2003")



## Panel models
library("plm")
pg <- plm.data(subset(Grunfeld, firm != "American Steel"), c("firm", "year"))

fm_fe <- plm(invest ~ value + capital, model = "within", data = pg)
summary(fm_fe)
coeftest(fm_fe, vcov = pvcovHC)

fm_reswar <- plm(invest ~ value + capital, data = pg,
  model = "random", random.method = "swar")
summary(fm_reswar)

## testing for random effects
fm_ols <- plm(invest ~ value + capital, data = pg, model = "pooling")
plmtest(fm_ols, type = "bp")
plmtest(fm_ols, type = "honda")

## Random effects models
fm_ream <- plm(invest ~ value + capital, data = pg, model = "random", random.method = "amemiya")
fm_rewh <- plm(invest ~ value + capital, data = pg, model = "random", random.method = "walhus")
fm_rener <- plm(invest ~ value + capital, data = pg, model = "random", random.method = "nerlove")

## Baltagi (2003), Tab. 2.1
rbind(
  "OLS(pooled)" = coef(fm_ols),
  "FE" = c(NA, coef(fm_fe)),
  "RE-SwAr" = coef(fm_reswar),
  "RE-Amemiya" = coef(fm_ream),
  "RE-WalHus" = coef(fm_rewh),
  "RE-Nerlove" = coef(fm_rener))

## Hausman test
phtest(fm_fe, fm_reswar)

## More examples can be found in:
## help("Baltagi2002")
## help("Greene2003")
}

\keyword{datasets}
