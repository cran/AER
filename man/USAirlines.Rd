\name{USAirlines}
\alias{USAirlines}

\title{Cost Data for US Airlines}

\description{
Cost data for six US airlines in 1970--1984.
}

\usage{data("USAirlines")}

\format{
A data frame containing 90 observations on 6 variables.
\describe{
  \item{firm}{factor indicating airline firm.}
  \item{year}{factor indicating year.}
  \item{output}{output revenue passenger miles index number.}
  \item{cost}{total cost (in USD 1000).}
  \item{price}{fuel price.}
  \item{load}{average capacity utilization of the fleet.}
}
}

\source{
Online complements to Greene (2003). Table F7.1.

\url{http://pages.stern.nyu.edu/~wgreene/Text/tables/tablelist5.htm}
}

\references{
Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Upper Saddle River, NJ: Prentice Hall.
}

\seealso{\code{\link{Greene2003}}}

\examples{
data("USAirlines")

## Example 7.2 in Greene (2003)
fm_full <- lm(log(cost) ~ log(output) + I(log(output)^2) + log(price) + load + year + firm,
  data = USAirlines)
fm_time <- lm(log(cost) ~ log(output) + I(log(output)^2) + log(price) + load + year,
  data = USAirlines)
fm_firm <- lm(log(cost) ~ log(output) + I(log(output)^2) + log(price) + load + firm,
  data = USAirlines)
fm_no <- lm(log(cost) ~ log(output) + I(log(output)^2) + log(price) + load, data = USAirlines)

## Table 7.2
anova(fm_full, fm_time)
anova(fm_full, fm_firm)
anova(fm_full, fm_no)

## alternatively, use plm()
library("plm")
usair <- plm.data(USAirlines, c("firm", "year"))
fm_full2 <- plm(log(cost) ~ log(output) + I(log(output)^2) + log(price) + load,
  data = usair, model = "within", effect = "twoways")
fm_time2 <- plm(log(cost) ~ log(output) + I(log(output)^2) + log(price) + load,
  data = usair, model = "within", effect = "time")
fm_firm2 <- plm(log(cost) ~ log(output) + I(log(output)^2) + log(price) + load,
  data = usair, model = "within", effect = "individual")
fm_no2 <- plm(log(cost) ~ log(output) + I(log(output)^2) + log(price) + load,
  data = usair, model = "pooling")
pFtest(fm_full2, fm_time2)
pFtest(fm_full2, fm_firm2)
pFtest(fm_full2, fm_no2)

## More examples can be found in:
## help("Greene2003")
}

\keyword{datasets}
