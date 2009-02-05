\name{Guns}
\alias{Guns}

\title{More Guns, Less Crime?}

\description{
Guns is a balanced panel of data on 50 US states, plus the District of Columbia (for a
total of 51 states), by year for 1977--1999.
}

\usage{data("Guns")}

\format{
A data frame containing 1,173 observations on 13 variables.
\describe{
  \item{state}{factor indicating state.}
  \item{year}{factor indicating year.}
  \item{violent}{violent crime rate (incidents per 100,000 members of the population).}
  \item{murder}{murder rate (incidents per 100,000).}
  \item{robbery}{robbery rate (incidents per 100,000).}
  \item{prisoners}{incarceration rate in the state in the previous year
    (sentenced prisoners per 100,000 residents; value for the previous year).}
  \item{afam}{percent of state population that is African-American, ages 10 to 64.}
  \item{cauc}{percent of state population that is Caucasian, ages 10 to 64.}
  \item{male}{percent of state population that is male, ages 10 to 29.}
  \item{population}{state population, in millions of people.}
  \item{income}{real per capita personal income in the state (US dollars).}
  \item{density}{population per square mile of land area, divided by 1,000.}
  \item{law}{factor. Does the state have a shall carry law in effect in that year?}
}
}

\details{
Each observation is a given state in a given year. 
There are a total of 51 states times 23 years = 1,173 observations.
}

\source{
Online complements to Stock and Watson (2007). 

\url{http://wps.aw.com/aw_stock_ie_2/0,12040,3332253-,00.html}
}

\references{
Ayres, I., and Donohue, J.J. (2003). Shooting Down the \sQuote{More Guns Less Crime} Hypothesis.
\emph{Stanford Law Review}, \bold{55}, 1193--1312.

Stock, J.H. and Watson, M.W. (2007). \emph{Introduction to Econometrics}, 2nd ed. Boston: Addison Wesley.
}

\seealso{\code{\link{StockWatson2007}}}

\examples{
## data
data("Guns")

## visualization
library("lattice")
xyplot(log(violent) ~ as.numeric(as.character(year)) | state, data = Guns, type = "l")

## Stock & Watson (2007), Empirical Exercise 10.1, pp. 376--377
fm1 <- lm(log(violent) ~ law, data = Guns)
coeftest(fm1, vcov = sandwich)

fm2 <- lm(log(violent) ~ law + prisoners + density + income + 
  population + afam + cauc + male, data = Guns)
coeftest(fm2, vcov = sandwich)

fm3 <- lm(log(violent) ~ law + prisoners + density + income + 
  population + afam + cauc + male + state, data = Guns)
printCoefmat(coeftest(fm3, vcov = sandwich)[1:9,])
            
fm4 <- lm(log(violent) ~ law + prisoners + density + income + 
  population + afam + cauc + male + state + year, data = Guns)
printCoefmat(coeftest(fm4, vcov = sandwich)[1:9,])
}

\keyword{datasets}
