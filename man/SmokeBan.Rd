\name{SmokeBan}
\alias{SmokeBan}

\title{Do Workplace Smoking Bans Reduce Smoking?}

\description{
Estimation of the effect of workplace smoking bans on smoking of indoor workers.
}

\usage{data("SmokeBan")}

\format{
A data frame containing 10,000 observations on 7 variables.
\describe{
  \item{smoker}{factor. Is the individual a current smoker?}
  \item{ban}{factor. Is there a work area smoking ban?}
  \item{age}{age in years.}
  \item{education}{factor indicating highest education level attained:
    high school (hs) drop out, high school graduate, some college,
    college graduate, master's degree (or higher).}
  \item{afam}{factor. Is the individual African-American?}  
  \item{hispanic}{factor. Is the individual Hispanic?} 
  \item{gender}{factor indicating gender.}
}
}

\details{
\code{SmokeBank} is a cross-sectional data set with observations on 10,000 indoor workers, which
is a subset of a 18,090-observation data set collected as part of the National Health
Interview Survey in 1991 and then again (with different respondents) in 1993.
The data set contains information on whether individuals were, or were not, subject to a workplace
smoking ban, whether or not the individuals smoked and other individual characteristics.
}

\source{
Online complements to Stock and Watson (2007). 

\url{http://wps.aw.com/aw_stock_ie_2/0,12040,3332253-,00.html}
}

\references{
Evans, W. N., Farrelly, M.C., and Montgomery, E. (1999). Do Workplace Smoking Bans Reduce Smoking?
\emph{American Economic Review}, \bold{89}, 728--747.

Stock, J.H. and Watson, M.W. (2007). \emph{Introduction to Econometrics}, 2nd ed. Boston: Addison Wesley.
}

\seealso{\code{\link{StockWatson2007}}}

\examples{
data("SmokeBan")

## proportion of non-smokers increases with education
plot(smoker ~ education, data = SmokeBan)

## proportion of non-smokers constant over age
plot(smoker ~ age, data = SmokeBan)
}

\keyword{datasets}
