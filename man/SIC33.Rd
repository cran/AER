\name{SIC33}
\alias{SIC33}
\title{SIC33 Production Data}
\description{
Statewide production data for primary metals industry (SIC 33).
}
\usage{data("SIC33")}

\format{
A data frame containing 27 observations on 3 variables.
\describe{
  \item{output}{Value added.}
  \item{labor}{Labor input.}
  \item{capital}{Capital stock.}
}
}

\source{
Online complements to Greene (2003). Table F6.1.

\url{https://pages.stern.nyu.edu/~wgreene/Text/tables/tablelist5.htm}
}

\references{
Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Upper Saddle River, NJ: Prentice Hall.
}

\seealso{\code{\link{Greene2003}}}

\examples{
\dontshow{ if(!requireNamespace("scatterplot3d")) {
  if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
    stop("not all packages required for the example are installed")
  } else q() }}
data("SIC33", package = "AER")

## Example 6.2 in Greene (2003)
## Translog model
fm_tl <- lm(output ~ labor + capital + I(0.5 * labor^2) + I(0.5 * capital^2) + I(labor * capital),
  data = log(SIC33))
## Cobb-Douglas model
fm_cb <- lm(output ~ labor + capital, data = log(SIC33))

## Table 6.2 in Greene (2003)
deviance(fm_tl)
deviance(fm_cb)
summary(fm_tl)
summary(fm_cb)
vcov(fm_tl)
vcov(fm_cb)

## Cobb-Douglas vs. Translog model
anova(fm_cb, fm_tl)
## hypothesis of constant returns
linearHypothesis(fm_cb, "labor + capital = 1")

## 3D Visualization
library("scatterplot3d")
  s3d <- scatterplot3d(log(SIC33)[,c(2, 3, 1)], pch = 16)
  s3d$plane3d(fm_cb, lty.box = "solid", col = 4)

## Interactive 3D Visualization
\donttest{
if(require("rgl")) {
  x <- log(SIC33)[,2]
  y <- log(SIC33)[,3]
  z <- log(SIC33)[,1]
  plot3d(x, y, z, type = "s", col = "gray", radius = 0.1)
  x <- seq(4.5, 7.5, by = 0.5)
  y <- seq(5.5, 10, by = 0.5)
  z <- outer(x, y, function(x, y) predict(fm_cb, data.frame(labor = x, capital = y)))
  surface3d(x, y, z, color = "blue", alpha = 0.5, shininess = 128)
}
}
}

\keyword{datasets}
