\name{TravelMode}
\alias{TravelMode}
\title{Travel Mode Choice Data}

\description{
Data on travel mode choice for travel between Sydney and Melbourne, Australia.
}

\usage{data("TravelMode")}

\format{
A data frame containing 840 observations on 4 modes for 210 individuals.
 \describe{
    \item{individual}{Factor indicating individual with levels \code{1} to \code{200}.}
    \item{mode}{Factor indicating travel mode with levels
      \code{"car"}, \code{"air"}, \code{"train"}, or \code{"bus"}.}
    \item{choice}{Factor indicating choice with levels \code{"no"} and \code{"yes"}.}
    \item{wait}{Terminal waiting time, 0 for car.}
    \item{vcost}{Vehicle cost component.}
    \item{travel}{Travel time in the vehicle.}
    \item{gcost}{Generalized cost measure.}
    \item{income}{Household income.}
    \item{size}{Party size.}
  }
}

\source{
Online complements to Greene (2003).

\url{http://pages.stern.nyu.edu/~wgreene/Text/tables/tablelist5.htm}
}

\references{
Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Upper Saddle River, NJ: Prentice Hall.
}

\seealso{\code{\link{Greene2003}}}

\examples{
data("TravelMode")

## overall proportions for chosen mode
with(TravelMode, prop.table(table(mode[choice == "yes"])))

## travel vs. waiting time for different travel modes
if(require("lattice")) {
  xyplot(travel ~ wait | mode, data = TravelMode)
}
}

\keyword{datasets}
