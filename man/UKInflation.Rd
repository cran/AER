\name{UKInflation}
\alias{UKInflation}
\title{UK Manufacturing Inflation Data}
\description{
Time series of observed and expected price changes in British manufacturing.
}
\usage{data("UKInflation")}
\format{
A quarterly multiple time series from 1972(1) to 1985(2) with 2 variables.
\describe{
\item{actual}{Actual inflation.}
\item{expected}{Expected inflation.}
}
}

\source{
Online complements to Greene (2003), Table F8.1.

\url{http://pages.stern.nyu.edu/~wgreene/Text/tables/tablelist5.htm}
}

\references{
Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Upper Saddle River, NJ: Prentice Hall.

Pesaran, M.H., and Hall, A.D. (1988). Tests of Non-nested Linear Regression Models
Subject To Linear Restrictions. \emph{Economics Letters}, \bold{27}, 341--348.
}

\seealso{\code{\link{Greene2003}}}

\examples{
data("UKInflation")
plot(UKInflation)
}

\keyword{datasets}
