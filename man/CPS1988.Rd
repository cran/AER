\name{CPS1988}
\alias{CPS1988}
\title{Determinants of Wages Data (CPS 1988)}

\description{
Cross-section data originating from the March 1988 Current Population Survey by the 
US Census Bureau.
}

\usage{data("CPS1988")}

\format{
A data frame containing 28,155 observations on 7 variables.
\describe{
  \item{wage}{Wage (in dollars per week).}
  \item{education}{Number of years of education.}
  \item{experience}{Number of years of potential work experience.}
  \item{ethnicity}{Factor with levels \code{"cauc"} and \code{"afam"} (African-American).}
  \item{smsa}{Factor. Does the individual reside in a Standard Metropolitan Statistical Area (SMSA)?}
  \item{region}{Factor with levels \code{"northeast"}, \code{"midwest"}, \code{"south"}, \code{"west"}.}
  \item{parttime}{Factor. Does the individual work part-time?}
}
}

\details{
A sample of men aged 18 to 70 with positive annual income greater than USD 50 in 1992, 
who are not self-employed nor working without pay. Wages are deflated by the deflator of 
Personal Consumption Expenditure for 1992.

A problem with CPS data is that it does not provide actual work experience. 
It is therefore customary to compute experience as \code{age - education - 6} 
(as was done by Bierens and Ginther, 2001), this may be considered potential experience. 
As a result, some respondents have negative experience.
}

\source{
\url{http://econ.la.psu.edu/~hbierens/MEDIAN.HTM}
}

\references{
Bierens, H.J., and Ginther, D. (2001). Integrated Conditional Moment Testing of Quantile Regression Models. 
\emph{Empirical Economics}, \bold{26}, 307--324.

Buchinsky, M. (1998). Recent Advances in Quantile Regression Models: A Practical Guide for Empirical Research. 
\emph{Journal of Human Resources}, \bold{33}, 88--126.
}

\seealso{\code{\link{CPS1985}}, \code{\link{CPSSW}}}

\examples{
## data and packages
library("quantreg")
data("CPS1988")
CPS1988$region <- relevel(CPS1988$region, ref = "south")

## Model equations: Mincer-type, quartic, Buchinsky-type
mincer <- log(wage) ~ ethnicity + education + experience + I(experience^2)
quart <- log(wage) ~ ethnicity + education + experience + I(experience^2) +
  I(experience^3) + I(experience^4)
buchinsky <- log(wage) ~ ethnicity * (education + experience + parttime) + 
  region*smsa + I(experience^2) + I(education^2) + I(education*experience)

## OLS and LAD fits (for LAD see Bierens and Ginter, Tables 1-3.A.)
mincer_ols <- lm(mincer, data = CPS1988)
mincer_lad <- rq(mincer, data = CPS1988)
quart_ols <- lm(quart, data = CPS1988)
quart_lad <- rq(quart, data = CPS1988)
buchinsky_ols <- lm(buchinsky, data = CPS1988)
buchinsky_lad <- rq(buchinsky, data = CPS1988)
}

\keyword{datasets}
