\name{Parade2005}
\alias{Parade2005}
\title{Parade Magazine 2005 Earnings Data}
\description{
US earnings data, as provided in an annual survey of Parade (here from 2005), the
Sunday newspaper magazine supplementing the Sunday (or Weekend) edition of many
daily newspapers in the USA. 
}
\usage{data("Parade2005")}

\format{
A data frame containing 130 observations on 5 variables.
\describe{
\item{earnings}{Annual personal earnings.}
\item{age}{Age in years.}
\item{gender}{Factor indicating gender.}
\item{state}{Factor indicating state.}
\item{celebrity}{Factor. Is the individual a celebrity?}
}
}

\details{
In addition to the four variables provided by Parade (earnings, age, gender, and state),
a fifth variable was introduced, the \dQuote{celebrity factor} (here actors, athletes, 
TV personalities, politicians, and CEOs are considered celebrities). The data are quite 
far from a simple random sample, there being substantial oversampling of celebrities.
}

\source{
Parade (2005). What People Earn. Issue March 13, 2005.
}

\examples{
## data
data("Parade2005")
attach(Parade2005)
summary(Parade2005)

## bivariate visualizations
plot(density(log(earnings), bw = "SJ"), type = "l", main = "log(earnings)")
rug(log(earnings))
plot(log(earnings) ~ gender, main = "log(earnings)")

## celebrity vs. non-celebrity earnings
noncel <- subset(Parade2005, celebrity == "no")
cel <- subset(Parade2005, celebrity == "yes")

library("ineq")
plot(Lc(noncel$earnings), main = "log(earnings)")
lines(Lc(cel$earnings), lty = 2)
lines(Lc(earnings), lty = 3)

Gini(noncel$earnings)
Gini(cel$earnings)
Gini(earnings)

## detach data
detach(Parade2005)
}

\keyword{datasets}
