\name{CPSSW}
\alias{CPSSW}
\alias{CPSSW9298}
\alias{CPSSW9204}
\alias{CPSSW04}
\alias{CPSSW3}
\alias{CPSSW8}
\alias{CPSSWEducation}

\title{Stock and Watson CPS Data Sets}

\description{
Stock and Watson (2007) provide several subsets created from 
March Current Population Surveys (CPS) with data on the
relationship of earnings and education over several year.
}

\usage{
data("CPSSW9204")
data("CPSSW9298")
data("CPSSW04")
data("CPSSW3")
data("CPSSW8")
data("CPSSWEducation")
}

\format{
\code{CPSSW9298}: A data frame containing 13,501 observations on 5 variables.
\code{CPSSW9204}: A data frame containing 15,588 observations on 5 variables.
\code{CPSSW04}: A data frame containing 7,986 observations on 4 variables.
\code{CPSSW3}: A data frame containing 20,999 observations on 3 variables.
\code{CPSSW8}: A data frame containing 61,395 observations on 5 variables.
\code{CPSSWEducation}: A data frame containing 2,950 observations on 4 variables.
\describe{
  \item{year}{factor indicating year.}
  \item{earnings}{average hourly earnings (sum of annual pretax wages, salaries, tips,
    and bonuses, divided by the number of hours worked annually).}
  \item{education}{number of years of education.}
  \item{degree}{factor indicating highest educational degree (\code{"bachelor"}
    or\code{"highschool"}).}
  \item{gender}{factor indicating gender.}
  \item{age}{age in years.} 
  \item{region}{factor indicating region of residence
   (\code{"Northeast"}, \code{"Midwest"}, \code{"South"}, \code{"West"}).}
}
}

\details{
Each month the Bureau of Labor Statistics in the US Department of Labor
conducts the Current Population Survey (CPS), which provides data
on labor force characteristics of the population, including the level of
employment, unemployment, and earnings. Approximately 65,000 randomly
selected US households are surveyed each month. The sample is chosen by
randomly selecting addresses from a database. 
Details can be found in the Handbook of Labor Statistics and is described
on the Bureau of Labor Statistics website (\url{http://www.bls.gov/}).

The survey conducted each March is more detailed than in other months
and asks questions about earnings during the previous year. 
The data sets contain data for 2004 (from the March 2005 survey), and some
also for earlier years (up to 1992). 

If education is given, it is for full-time workers, defined as
workers employed more than 35 hours per week for at least 48 weeks in the
previous year. Data are provided for workers whose 
highest educational achievement is a high school diploma and a bachelor's degree.

Earnings for years earlier than 2004 were adjusted for inflation by putting them
in 2004 USD using the Consumer Price Index (CPI). From 1992 to 2004, the
price of the CPI market basket rose by 34.6\%. To make earnings in 1992 and 2004
comparable, 1992 earnings are inflated by the amount of overall CPI price
inflation, by multiplying 1992 earnings by 1.346 to put them into 2004 dollars.

\code{CPSSW9204} provides the distribution of earnings in the US in 1992 and 2004
for college-educated full-time workers aged 25--34.
\code{CPSSW04} is a subset of \code{CPSSW9204} and provides the distribution of
earnings in the US in 2004 for college-educated full-time workers aged 25--34.
\code{CPSSWEducation} is similar (but not a true subset) and contains the
distribution of earnings in the US in 2004 for college-educated full-time workers
aged 29--30.
\code{CPSSW8} contains a larger sample with workers aged 21--64, additionally
providing information about the region of residence.
\code{CPSSW9298} is similar to \code{CPSSW9204} providing data from 1992 and 1998
(with the 1992 subsets not being exactly identical).
\code{CPSSW3} provides trends (from 1992 to 2004) in hourly earnings in the US of
working college graduates aged 25--34 (in 2004 USD).
}

\source{
Online complements to Stock and Watson (2007). 

\url{http://wps.aw.com/aw_stock_ie_2/}
}

\references{
Stock, J.H. and Watson, M.W. (2007). \emph{Introduction to Econometrics}, 2nd ed. Boston: Addison Wesley.
}

\seealso{\code{\link{StockWatson2007}}, \code{\link{CPS1985}}, \code{\link{CPS1988}}}

\examples{
data("CPSSW3")
with(CPSSW3, interaction.plot(year, gender, earnings))

## Stock and Watson, p. 165
data("CPSSWEducation")
plot(earnings ~ education, data = CPSSWEducation)
fm <- lm(earnings ~ education, data = CPSSWEducation)
coeftest(fm, vcov = sandwich)
abline(fm)
}

\keyword{datasets}
