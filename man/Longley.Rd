\name{Longley}
\alias{Longley}
\title{Longley's Regression Data}
\description{
US macroeconomic time series, 1947--1962.
}
\usage{data("Longley")}
\format{
An annual multiple time series from 1947 to 1962 with 4 variables.
 \describe{
    \item{employment}{Number of people employed (in 1000s).}
    \item{price}{GNP deflator.}
    \item{gnp}{Gross national product.}
    \item{armedforces}{Number of people in the armed forces.}
 }
}

\details{
An extended version of this data set, formatted as a \code{"data.frame"}
is available as \code{\link[datasets]{longley}} in base R.
}

\source{
Online complements to Greene (2003). Table F4.2.

\url{http://pages.stern.nyu.edu/~wgreene/Text/tables/tablelist5.htm}
}

\references{
Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Upper Saddle River, NJ: Prentice Hall.

Longley, J.W. (1967). An Appraisal of Least-Squares Programs from the Point of View of the User.
\emph{Journal of the American Statistical Association}, \bold{62}, 819--841.
}

\seealso{\code{\link[datasets]{longley}}, \code{\link{Greene2003}}}

\examples{
data("Longley")
library("dynlm")

## Example 4.6 in Greene (2003)
fm1 <- dynlm(employment ~ time(employment) + price + gnp + armedforces,
  data = Longley)
fm2 <- update(fm1, end = 1961)
cbind(coef(fm2), coef(fm1))

## Figure 4.3 in Greene (2003)
plot(rstandard(fm2), type = "b", ylim = c(-3, 3))
abline(h = c(-2, 2), lty = 2)
}

\keyword{datasets}
