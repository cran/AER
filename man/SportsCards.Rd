\name{SportsCards}
\alias{SportsCards}

\title{Endowment Effect for Sports Cards}

\description{
Trading sports cards: Does ownership increase the value of goods to consumers?
}

\usage{data("SportsCards")}

\format{
A data frame containing 148 observations on 9 variables.
\describe{
  \item{good}{factor. Was the individual given good A or B (see below)?}
  \item{dealer}{factor. Was the individual a dealer?}
  \item{permonth}{number of trades per month reported by the individual.}
  \item{years}{number of years that the individual has been trading.}
  \item{income}{factor indicating income group (in 1000 USD).}
  \item{gender}{factor indicating gender.}
  \item{education}{factor indicating highest level of education (8th grade or less,
    high school, 2-year college, other post-high school, 4-year college or graduate school).}
  \item{age}{age in years.}
  \item{trade}{factor. Did the individual trade the good he was given for the other good?}
}
}

\details{
\code{SportsCards} contains data from 148 randomly selected traders who attended
a trading card show in Orlando, Florida, in 1998. Traders were randomly given one
of two sports collectables, say good A or good B, that had approximately equal market
value. Those receiving good A were then given the option of trading good A for good B
with the experimenter; those receiving good B were given the option of trading good B
for good A with the experimenter. Good A was a ticket stub from the game that Cal Ripken Jr.
set the record for consecutive games played, and Good B was a souvenir
from the game that Nolan Ryan won his 300th game.
}

\source{
Online complements to Stock and Watson (2007). 

\url{http://wps.aw.com/aw_stock_ie_2/}
}

\references{
List, J.A. (2003). Does Market Experience Eliminate Market Anomalies?
\emph{Quarterly Journal of Economcis}, \bold{118}, 41--71.

Stock, J.H. and Watson, M.W. (2007). \emph{Introduction to Econometrics}, 2nd ed. Boston: Addison Wesley.
}

\seealso{\code{\link{StockWatson2007}}}

\examples{
data("SportsCards")
summary(SportsCards)

plot(trade ~ permonth, data = SportsCards,
  ylevels = 2:1, breaks = c(0, 5, 10, 20, 30, 70))
plot(trade ~ years, data = SportsCards,
  ylevels = 2:1, breaks = c(0, 5, 10, 20, 60))
}

\keyword{datasets}
