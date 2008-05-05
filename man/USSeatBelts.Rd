\name{USSeatBelts}
\alias{USSeatBelts}

\title{Effects of Mandatory Seat Belt Laws in the US}

\description{
Balanced panel data for the years 1983--1997 from 50 US States, plus the District of Columbia,
for assessing traffic fatalities and seat belt usage.
}

\usage{data("USSeatBelts")}

\format{
A data frame containing 765 observations on 12 variables.
\describe{
  \item{state}{factor indicating US state (abbreviation).}
  \item{year}{factor indicating year.}
  \item{miles}{millions of traffic miles per year.}
  \item{fatalities}{number of fatalities per million of traffic miles
    (absolute frequencies of fatalities = \code{fatalities} times \code{miles}).}
  \item{seatbelt}{seat belt usage rate, as self-reported by state population surveyed.}  
  \item{speed65}{factor. Is there a 65 mile per hour speed limit?}
  \item{speed70}{factor. Is there a 70 (or higher) mile per hour speed limit?}  
  \item{drinkage}{factor. Is there a minimum drinking age of 21 years?} 
  \item{alcohol}{factor. Is there a maximum of 0.08 blood alcohol content?}
  \item{income}{median per capita income (in current US dollar).}
  \item{age}{mean age.}
  \item{enforce}{factor indicating seat belt law enforcement
    (\code{"no"}, \code{"primary"}, \code{"secondary"}).}
}
}

\details{
Some data series from Cohen and Einav (2003) have not been included in the data frame.
}

\source{
Online complements to Stock and Watson (2007). 

\url{http://wps.aw.com/aw_stock_ie_2/}
}

\references{
Cohen, A., and Einav, L. (2003). The Effects of Mandatory Seat Belt Laws on Driving
Behavior and Traffic Fatalities.
\emph(The Review of Economics and Statistics), \bold(85), 828--843

Stock, J.H. and Watson, M.W. (2007). \emph{Introduction to Econometrics}, 2nd ed. Boston: Addison Wesley.
}

\seealso{\code{\link{StockWatson2007}}}

\examples{
data("USSeatBelts")
summary(USSeatBelts)

library("lattice")
xyplot(fatalities ~ as.numeric(as.character(year)) | state, data = USSeatBelts, type = "l")
}

\keyword{datasets}
