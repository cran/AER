\name{EuroEnergy}
\alias{EuroEnergy}
\title{European Energy Consumption Data}
\description{
Cross-section data on energy consumption for 20 European countries, for the year 1980. 
}
\usage{data("EuroEnergy")}
\format{
A data frame containing 20 observations on 2 variables. 
\describe{
\item{gdp}{Real gross domestic product for the year 1980 (in million 1975 US dollars).}
\item{energy}{Aggregate energy consumption (in million kilograms coal equivalence).}
}
}

\source{
The data are from Baltagi (2002) and available at

\url{http://www.springeronline.com/sgw/cda/frontpage/0,10735,4-165-2-107420-0,00.html}
}

\references{
Baltagi, B.H. (2002). \emph{Econometrics}, 3rd ed. Berlin, Springer. 
}

\seealso{\code{\link{Baltagi2002}}}

\examples{
data("EuroEnergy")
energy_lm <- lm(log(energy) ~ log(gdp), data = EuroEnergy)
influence.measures(energy_lm)
}

\keyword{datasets}
