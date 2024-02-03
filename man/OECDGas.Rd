\name{OECDGas}
\alias{OECDGas}
\title{Gasoline Consumption Data}
\description{
Panel data on gasoline consumption in 18 OECD countries over 19 years, 1960--1978. 
}
\usage{data("OECDGas")}
\format{
A data frame containing 342 observations on 6 variables.
\describe{
  \item{country}{Factor indicating country.}
  \item{year}{Year.}
  \item{gas}{Logarithm of motor gasoline consumption per car.}
  \item{income}{Logarithm of real per-capita income.}
  \item{price}{Logarithm of real motor gasoline price.}
  \item{cars}{Logarithm of the stock of cars per-capita.}
}
}

\source{
The data is from Baltagi (2002).
}

\references{
Baltagi, B.H. (2002). \emph{Econometrics}, 3rd ed. Berlin, Springer. 

Baltagi, B.H. and Griffin, J.M. (1983). Gasoline Demand in the OECD: An Application
of Pooling and Testing Procedures. \emph{European Economic Review}, \bold{22}, 117--137.
}

\seealso{\code{\link{Baltagi2002}}}

\examples{
\dontshow{ if(!requireNamespace("lattice")) {
  if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
    stop("not all packages required for the example are installed")
  } else q() }}
data("OECDGas")

library("lattice")
xyplot(exp(cars) ~ year | country, data = OECDGas, type = "l")
xyplot(exp(gas) ~ year | country, data = OECDGas, type = "l")
}

\keyword{datasets}
