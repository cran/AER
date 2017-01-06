\name{MotorCycles2}
\alias{MotorCycles2}
\title{Motor Cycles in The Netherlands}
\description{
Time series of stock of motor cycles (two wheels) in The Netherlands
(in thousands).
}
\usage{data("MotorCycles2")}

\format{
An annual univariate time series from 1946 to 2012.
}

\details{This is an update of the series that was available with Franses (1998). However, the values for the years 1992 and 1993 differ.}

\source{
Online complements to Franses, van Dijk and Opschoor (2014).

\url{http://www.cambridge.org/us/academic/subjects/economics/econometrics-statistics-and-mathematical-economics/time-series-models-business-and-economic-forecasting-2nd-edition}
}

\references{
Franses, P.H. (1998). \emph{Time Series Models for Business and Economic Forecasting}.
Cambridge, UK: Cambridge University Press.

Franses, P.H., van Dijk, D. and Opschoor, A. (2014). \emph{Time Series Models for Business and Economic Forecasting}, 2nd ed.
Cambridge, UK: Cambridge University Press.
}

\seealso{\code{\link{Franses1998}}, \code{\link{MotorCycles}}}

\examples{
data("MotorCycles2")
plot(MotorCycles2)
}

\keyword{datasets}
