\name{ManufactCosts}
\alias{ManufactCosts}
\title{Manufacturing Costs Data}
\description{
US time series data on prices and cost shares in manufacturing, 1947--1971.
}
\usage{data("ManufactCosts")}
\format{
An annual multiple time series from 1947 to 1971 with 9 variables.
\describe{
  \item{cost}{Cost index.}
  \item{capitalcost}{Capital cost share.}
  \item{laborcost}{Labor cost share.}
  \item{energycost}{Energy cost share.}
  \item{materialscost}{Materials cost share.}
  \item{capitalprice}{Capital price.}
  \item{laborprice}{Labor price.}
  \item{energyprice}{Energy price.}
  \item{materialsprice}{Materials price.}
}
}

\source{
Online complements to Greene (2003).

\url{http://pages.stern.nyu.edu/~wgreene/Text/tables/tablelist5.htm}
}

\references{
Berndt, E. and Wood, D. (1975). Technology, Prices, and the Derived Demand for Energy.
\emph{Review od Economics and Statistics}, \bold{57}, 376--384.

Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Upper Saddle River, NJ: Prentice Hall.
}

\seealso{\code{\link{Greene2003}}}

\examples{
data("ManufactCosts")
plot(ManufactCosts)
}

\keyword{datasets}
