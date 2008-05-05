\name{WeakInstrument}
\alias{WeakInstrument}

\title{Artificial Weak Instrument Data}

\description{
Artificial data set to illustrate the problem of weak instruments.
}

\usage{data("WeakInstrument")}

\format{
A data frame containing 200 observations on 3 variables.
\describe{
  \item{y}{dependent variable.}
  \item{x}{regressor variable.}
  \item{z}{instrument variable.}
}
}

\source{
Online complements to Stock and Watson (2007). 

\url{http://wps.aw.com/aw_stock_ie_2/0,12040,3332253-,00.html}
}

\references{
Stock, J.H. and Watson, M.W. (2007). \emph{Introduction to Econometrics}, 2nd ed. Boston: Addison Wesley.
}

\seealso{\code{\link{StockWatson2007}}}

\examples{
data("WeakInstrument")
fm <- ivreg(y ~ x | z, data = WeakInstrument)
summary(fm)
}

\keyword{datasets}
