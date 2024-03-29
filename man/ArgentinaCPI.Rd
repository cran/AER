\name{ArgentinaCPI}
\alias{ArgentinaCPI}
\title{Consumer Price Index in Argentina}
\description{
Time series of consumer price index (CPI) in Argentina (index with 1969(4) = 1).
}
\usage{data("ArgentinaCPI")}

\format{
A quarterly univariate time series from 1970(1) to 1989(4).
}

\source{
Online complements to Franses (1998).
}

\references{
De Ruyter van Steveninck, M.A. (1996). \emph{The Impact of Capital Imports; Argentina 1970--1989}.
Amsterdam: Thesis Publishers.

Franses, P.H. (1998). \emph{Time Series Models for Business and Economic Forecasting}.
Cambridge, UK: Cambridge University Press.
}

\seealso{\code{\link{Franses1998}}}

\examples{
\dontshow{ if(!requireNamespace("dynlm")) {
  if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
    stop("not all packages required for the example are installed")
  } else q() }}
data("ArgentinaCPI")
plot(ArgentinaCPI)
plot(log(ArgentinaCPI))

library("dynlm")
## estimation sample 1970.3-1988.4 means
acpi <- window(ArgentinaCPI, start = c(1970,1), end = c(1988,4)) 

## eq. (3.90), p.54
acpi_ols <- dynlm(d(log(acpi)) ~ L(d(log(acpi))))
summary(acpi_ols)

## alternatively
ar(diff(log(acpi)), order.max = 1, method = "ols") 
}

\keyword{datasets}
