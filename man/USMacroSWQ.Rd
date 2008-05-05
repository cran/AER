\name{USMacroSWQ}
\alias{USMacroSWQ}

\title{Quarterly US Macroeconomic Data (1947--2004, Stock \& Watson)}

\description{
Time-series data on 2 US macroeconomic variables for 1947--2004.
}
\usage{data("USMacroSWQ")}

\format{
A quarterly multiple time series from 1947(1) to 2004(4) with 2 variables.
\describe{
  \item{gdp}{real GDP for the United States in billions of chained (2000) dollars
    seasonally adjusted, annual rate.}
  \item{tbill}{3-month treasury bill rate. Quarterly averages of daily dates
    in percentage points at an annual rate.}
  }
}

\source{
Online complements to Stock and Watson (2007). 

\url{http://wps.aw.com/aw_stock_ie_2/}
}

\references{
Stock, J.H. and Watson, M.W. (2007). \emph{Introduction to Econometrics}, 2nd ed. Boston: Addison Wesley.
}

\seealso{\code{\link{StockWatson2007}}, \code{\link{USMacroSW}}, \code{\link{USMacroSWM}},
  \code{\link{USMacroB}}, \code{\link{USMacroG}}}

\examples{
data("USMacroSWQ")
plot(USMacroSWQ)
}

\keyword{datasets}
