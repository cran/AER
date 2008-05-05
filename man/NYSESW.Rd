\name{NYSESW}
\alias{NYSESW}

\title{Daily NYSE Composite Index}

\description{
A daily time series from 1990 to 2005 of the New York Stock Exchange composite index.
}

\usage{data("NYSESW")}

\format{
A daily univariate time series from 1990-01-02 to 2005-11-11 (of class
\code{"zoo"} with \code{"Date"} index).
}

\source{
Online complements to Stock and Watson (2007). 

\url{http://wps.aw.com/aw_stock_ie_2/}
}

\references{
Stock, J.H. and Watson, M.W. (2007). \emph{Introduction to Econometrics}, 2nd ed. Boston: Addison Wesley.
}

\seealso{\code{\link{StockWatson2007}}}

\examples{
## returns
data("NYSESW")
ret <- 100 * diff(log(NYSESW))
plot(ret)

## Stock and Watson (2007), p. 667, GARCH(1,1) model
library("tseries")
fm <- garch(coredata(ret))
summary(fm)
}


\keyword{datasets}
