\name{Fertility}
\alias{Fertility}
\alias{Fertility2}

\title{Fertility and Women's Labor Supply}

\description{
Cross-section data from the 1980 US Census on married women
aged 21--35 with two or more children.
}

\usage{
data("Fertility")
data("Fertility2")
}

\format{
A data frame containing 254,654 (and 30,000, respectively) observations on 8 variables.
\describe{
  \item{morekids}{factor. Does the mother have more than 2 children?}
  \item{gender1}{factor indicating gender of first child.}
  \item{gender2}{factor indicating gender of second child.}
  \item{age}{age of mother at census.}
  \item{afam}{factor. Is the mother African-American?}  
  \item{hispanic}{factor. Is the mother Hispanic?} 
  \item{other}{factor. Is the mother's ethnicity neither African-American nor
    Hispanic, nor Caucasian? (see below)}
  \item{work}{number of weeks in which the mother worked in 1979.}
}
}

\details{
\code{Fertility2} is a random subset of \code{Fertility} with 30,000 observations.

There are conflicts in the ethnicity coding (see also examples). Hence, it was
not possible to create a single factor and the original three indicator
variables have been retained.

Not all variables from Angrist and Evans (1998) have been included.
}

\source{
Online complements to Stock and Watson (2007). 

\url{http://wps.aw.com/aw_stock_ie_2/0,12040,3332253-,00.html}
}

\references{
Angrist, J.D., and Evans, W.N. (1998). Children and Their Parents' Labor Supply: Evidence from Exogenous Variation in Family Size
\emph{American Economic Review}, \bold{88}, 450--477.

Stock, J.H. and Watson, M.W. (2007). \emph{Introduction to Econometrics}, 2nd ed. Boston: Addison Wesley.
}

\seealso{\code{\link{StockWatson2007}}}

\examples{
data("Fertility2")

## conflicts in ethnicity coding
ftable(xtabs(~ afam + hispanic + other, data = Fertility2))

## create convenience variables
Fertility2$mkids <- with(Fertility2, as.numeric(morekids) - 1)
Fertility2$samegender <- with(Fertility2, factor(gender1 == gender2))
Fertility2$twoboys <- with(Fertility2, factor(gender1 == "male" & gender2 == "male"))
Fertility2$twogirls <- with(Fertility2, factor(gender1 == "female" & gender2 == "female"))

## similar to Angrist and Evans, p. 462
fm1 <- lm(mkids ~ samegender, data = Fertility2)
summary(fm1)

fm2 <- lm(mkids ~ gender1 + gender2 + samegender + age + afam + hispanic + other, data = Fertility2)
summary(fm2)

fm3 <- lm(mkids ~ gender1 + twoboys + twogirls + age + afam + hispanic + other, data = Fertility2)
summary(fm3)
}

\keyword{datasets}
