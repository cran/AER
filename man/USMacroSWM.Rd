\name{USMacroSWM}
\alias{USMacroSWM}

\title{Monthly US Macroeconomic Data (1947--2004, Stock \& Watson)}

\description{
Time-series data on 4 US macroeconomic variables for 1947--2004.
}
\usage{data("USMacroSWM")}

\format{
A monthly multiple time series from 1947(1) to 2004(4) with 4 variables.
\describe{
  \item{production}{index of industrial production.}
  \item{oil}{oil price shocks, starting 1948(1).}
  \item{cpi}{all-items consumer price index.}
  \item{expenditure}{personal consumption expenditures price deflator, starting 1959(1).}
  }
}

\source{
Online complements to Stock and Watson (2007). 

\url{http://wps.aw.com/aw_stock_ie_2/0,12040,3332253-,00.html}
}

\references{
Stock, J.H. and Watson, M.W. (2007). \emph{Introduction to Econometrics}, 2nd ed. Boston: Addison Wesley.
}

\seealso{\code{\link{StockWatson2007}}, \code{\link{USMacroSW}}, \code{\link{USMacroSWQ}},
  \code{\link{USMacroB}}, \code{\link{USMacroG}}}

\examples{
data("USMacroSWM")
plot(USMacroSWM)
}

\keyword{datasets}
