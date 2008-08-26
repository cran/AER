\name{USMoney}
\alias{USMoney}

\title{USMoney}

\description{
Money, output and price deflator time series data, 1950--1983.  
}

\usage{data("USMoney")}

\format{
A quarterly multiple time series from 1950 to 1983 with 3 variables.
\describe{
  \item{gnp}{nominal GNP.}
  \item{m1}{M1 measure of money stock.}
  \item{deflator}{implicit price deflator for GNP.}
}
}

\source{
Online complements to Greene (2003), Table F20.2.

\url{http://pages.stern.nyu.edu/~wgreene/Text/tables/tablelist5.htm}
}

\references{
Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Upper Saddle River, NJ: Prentice Hall.
}

\seealso{\code{\link{Greene2003}}}

\examples{
data("USMoney")
plot(USMoney)
}

\keyword{datasets}
