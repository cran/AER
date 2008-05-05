\name{HousePrices}
\alias{HousePrices}

\title{House Prices in the City of Windsor, Canada}

\description{
Sales prices of houses sold in the city of Windsor, Canada, during July, August and September, 1987.
}

\usage{data("HousePrices")}

\format{
A data frame containing 546 observations on 12 variables.
\describe{
  \item{price}{Sale price of a house.}
  \item{lotsize}{Lot size of a property in square feet.}
  \item{bedrooms}{Number of bedrooms.}
  \item{bathrooms}{Number of full bathrooms.}
  \item{stories}{Number of stories excluding basement.}
  \item{driveway}{Factor. Does the house have a driveway?}
  \item{recreation}{Factor. Does the house have a recreational room?}
  \item{fullbase}{Factor. Does the house have a full finished basement?}
  \item{gasheat}{Factor. Does the house use gas for hot water heating?}
  \item{aircon}{Factor. Is there central air conditioning?}
  \item{garage}{Number of garage places.}
  \item{prefer}{Factor. Is the house located in the preferred neighbourhood of the city?}
}
}

\source{
Journal of Applied Econometrics Data Archive. 

\url{http://www.econ.queensu.ca/jae/1996-v11.6/anglin-gencay/}
}

\references{
Anglin, P., and Gencay, R. (1996). Semiparametric Estimation of a Hedonic Price Function. 
\emph{Journal of Applied Econometrics}, \bold{11}, 633--648.

Verbeek, M. (2004). \emph{A Guide to Modern Econometrics}, 2nd ed. Chichester, UK: John Wiley.
}

\examples{
data("HousePrices")

### Anglin + Gencay (1996), Table II
fm_ag <- lm(log(price) ~ driveway + recreation + fullbase + gasheat + 
  aircon + garage + prefer + log(lotsize) + log(bedrooms) + 
  log(bathrooms) + log(stories), data = HousePrices)

### Anglin + Gencay (1996), Table III
fm_ag2 <- lm(log(price) ~ driveway + recreation + fullbase + gasheat + 
  aircon + garage + prefer + log(lotsize) + bedrooms + 
  bathrooms + stories, data = HousePrices)

### Verbeek (2004), Table 3.1
fm <- lm(log(price) ~ log(lotsize) + bedrooms + bathrooms + aircon, data = HousePrices)
summary(fm)

### Verbeek (2004), Table 3.2
fm_ext <- lm(log(price) ~ . - lotsize + log(lotsize), data = HousePrices)
summary(fm_ext)

### Verbeek (2004), Table 3.3
fm_lin <- lm(price ~ . , data = HousePrices)
summary(fm_lin)
}

\keyword{datasets}
