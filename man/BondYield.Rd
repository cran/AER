\name{BondYield}
\alias{BondYield}
\title{Bond Yield Data}

\description{
Monthly averages of the yield on a Moody's Aaa rated corporate bond (in percent/year).  
}

\usage{data("BondYield")}
\format{
A monthly univariate time series from 1990(1) to 1994(12).
}

\source{
Online complements to Greene (2003), Table F20.1.

\url{http://pages.stern.nyu.edu/~wgreene/Text/tables/tablelist5.htm}
}

\references{
Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Upper Saddle River, NJ: Prentice Hall.
}

\seealso{\code{\link{Greene2003}}}

\examples{
data("BondYield")
plot(BondYield)
}

\keyword{datasets}
