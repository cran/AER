\name{GrowthSW}
\alias{GrowthSW}

\title{Determinants of Economic Growth}

\description{
Data on average growth rates over 1960--1995 for 65 countries, along with
variables that are potentially related to growth.
}

\usage{data("GrowthSW")}

\format{
A data frame containing 65 observations on 6 variables.
\describe{
  \item{growth}{average annual percentage growth of real GDP from 1960 to 1995.}
  \item{rgdp60}{value of GDP per capita in 1960, converted to 1960 US dollars.}
  \item{tradeshare}{average share of trade in the economy from 1960 to 1995, 
        measured as the sum of exports (X) plus imports (M), divided by GDP;
        that is, the average value of (X + M)/GDP from 1960 to 1995.}
  \item{education}{average number of years of schooling of adult residents in that country in 1960.}
  \item{revolutions}{average annual number of revolutions, 
        insurrections (successful or not) and coup d'etats in that country from 1960 to 1995.}
  \item{assassinations}{average annual number of political assassinations 
        in that country from 1960 to 1995 (in per million population).}
}
}

\source{
Online complements to Stock and Watson (2007). 

\url{http://wps.aw.com/aw_stock_ie_2}
}

\references{
Beck, T., Levine, R., and Loayza, N. (2000). Finance and the Sources of Growth.
\emph{Journal of Financial Economics}, \bold{58}, 261--300.

Stock, J. H. and Watson, M. W. (2007). \emph{Introduction to Econometrics}, 2nd ed. Boston: Addison Wesley.
}

\seealso{\code{\link{StockWatson2007}}, \code{\link{GrowthDJ}}, \code{\link{OECDGrowth}}}

\examples{
data("GrowthSW")
summary(GrowthSW)
}

\keyword{datasets}
