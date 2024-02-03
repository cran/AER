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
    \item{individual}{Factor indicating individual with levels \code{1} to \code{210}.}
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

\url{https://pages.stern.nyu.edu/~wgreene/Text/tables/tablelist5.htm}
}

\references{
Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Upper Saddle River, NJ: Prentice Hall.
}

\seealso{\code{\link{Greene2003}}}

\examples{
\dontshow{ if(!requireNamespace("lattice") ||
              !requireNamespace("mlogit")) {
  if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
    stop("not all packages required for the example are installed")
  } else q() }}
data("TravelMode", package = "AER")

## overall proportions for chosen mode
with(TravelMode, prop.table(table(mode[choice == "yes"])))

## travel vs. waiting time for different travel modes
library("lattice")
xyplot(travel ~ wait | mode, data = TravelMode)

## Greene (2003), Table 21.11, conditional logit model
library("mlogit")
TravelMode$incair <- with(TravelMode, income * (mode == "air"))
tm_cl <- mlogit(choice ~ gcost + wait + incair, data = TravelMode,
  shape = "long", alt.var = "mode", reflevel = "car")
summary(tm_cl)
}

\keyword{datasets}
