\name{DJIA8012}
\alias{DJIA8012}
\title{Dow Jones Industrial Average (DJIA) index}

\description{
Time series of the Dow Jones Industrial Average (DJIA) index.
}

\usage{data("DJIA8012")}

\format{
A daily univariate time series from 1980-01-01 to 2012-12-31 (of class \code{"zoo"} with \code{"Date"} index).
}

\source{
Online complements to Franses, van Dijk and Opschoor (2014).

\url{https://www.cambridge.org/us/academic/subjects/economics/econometrics-statistics-and-mathematical-economics/time-series-models-business-and-economic-forecasting-2nd-edition}
}

\references{
Franses, P.H., van Dijk, D. and Opschoor, A. (2014). \emph{Time Series Models for Business and Economic Forecasting}, 2nd ed.
Cambridge, UK: Cambridge University Press.
}

\examples{
data("DJIA8012")
plot(DJIA8012)

# p.26, Figure 2.18
dldjia <- diff(log(DJIA8012))
plot(dldjia)

# p.141, Figure 6.4
plot(window(dldjia, start = "1987-09-01", end = "1987-12-31"))

# p.167, Figure 7.1
dldjia9005 <- window(dldjia, start = "1990-01-01", end = "2005-12-31")
qqnorm(dldjia9005)
qqline(dldjia9005, lty = 2)

# p.170, Figure 7.4
acf(dldjia9005,  na.action = na.exclude, lag.max = 250, ylim =  c(-0.1, 0.25))
acf(dldjia9005^2,  na.action = na.exclude, lag.max = 250, ylim =  c(-0.1, 0.25))
acf(abs(dldjia9005),  na.action = na.exclude, lag.max = 250, ylim =  c(-0.1, 0.25))
}

\keyword{datasets}
