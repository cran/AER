\name{MotorCycles}
\alias{MotorCycles}
\title{Motor Cycles in The Netherlands}
\description{
Time series of stock of motor cycles (two wheels) in The Netherlands
(in thousands).
}
\usage{data("MotorCycles")}

\format{
An annual univariate time series from 1946 to 1993.
}

\details{An updated version is available under the name \code{MotorCycles2}. However, the values for the years 1992 and 1993 differ there.}

\source{
Online complements to Franses (1998).
}

\references{
Franses, P.H. (1998). \emph{Time Series Models for Business and Economic Forecasting}.
Cambridge, UK: Cambridge University Press.
}

\seealso{\code{\link{Franses1998}}, \code{\link{MotorCycles2}}}

\examples{
data("MotorCycles")
plot(MotorCycles)
}

\keyword{datasets}
