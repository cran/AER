\name{ShipAccidents}
\alias{ShipAccidents}

\title{Ship Accidents}

\description{
Data on ship accidents.
}

\usage{data("ShipAccidents")}

\format{
A data frame containing 40 observations on 5 ship types in 4 vintages and 2 service periods. 
 \describe{
    \item{type}{factor with levels \code{"A"} to \code{"E"} for the different ship types,}
    \item{construction}{factor with levels \code{"1960-64"}, \code{"1965-69"}, \code{"1970-74"}, 
      \code{"1975-79"} for the periods of construction,}
    \item{operation}{factor with levels \code{"1960-74"}, \code{"1975-79"} for the periods of
      operation,}
    \item{service}{aggregate months of service,}
    \item{incidents}{number of damage incidents.}
    }
}

\details{
The data are from McCullagh and Nelder (1989, p. 205, Table 6.2) and were
also used by Greene (2003, Ch. 21), see below.

There are five ships (observations 7, 15, 23, 31, 39) with an operation period
\emph{before} the construction period, hence the variables \code{service} and
\code{incidents} are necessarily 0. An additional observation (34) has entries
representing \emph{accidentally empty cells} (see McCullagh and Nelder, 1989, p. 205).

It is a bit unclear what exactly the above means. In any case, the models are fit
only to those observations with \code{service > 0}.
}

\source{
Online complements to Greene (2003).

\url{http://pages.stern.nyu.edu/~wgreene/Text/tables/tablelist5.htm}
}

\references{
Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Upper Saddle River, NJ: Prentice Hall.

McCullagh, P. and Nelder, J.A. (1989). \emph{Generalized Linear Models}, 2nd edition.
London: Chapman \& Hall.
}

\seealso{\code{\link{Greene2003}}}

\examples{
data("ShipAccidents")
sa <- subset(ShipAccidents, service > 0)

## Greene (2003), Table 21.20
## (see also McCullagh and Nelder, 1989, Table 6.3)
sa_full <- glm(incidents ~ type + construction + operation, family = poisson,
  data = sa, offset = log(service))
summary(sa_full)

sa_notype <- glm(incidents ~ construction + operation, family = poisson,
  data = sa, offset = log(service))
summary(sa_notype)

sa_noperiod <- glm(incidents ~ type + operation, family = poisson,
  data = sa, offset = log(service))
summary(sa_noperiod)

## model comparison
anova(sa_full, sa_notype, test = "Chisq")
anova(sa_full, sa_noperiod, test = "Chisq")

## test for overdispersion
dispersiontest(sa_full)
dispersiontest(sa_full, trafo = 2)
}

\keyword{datasets}
