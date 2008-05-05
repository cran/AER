\name{HealthInsurance}
\alias{HealthInsurance}

\title{Medical Expenditure Panel Survey Data}

\description{
Cross-section data originating from the Medical Expenditure Panel Survey survey conducted in 1996.
}

\usage{data("HealthInsurance")}

\format{
A data frame containing 8,802 observations on 11 variables.
\describe{
  \item{health}{factor. Is the self-reported health status \dQuote{healthy}?.}
  \item{age}{age in years.}
  \item{limit}{factor. Is there any limitation?}
  \item{gender}{factor indicating gender.}
  \item{insurance}{factor. Does the individual have a health insurance?}
  \item{married}{factor. Is the individual married?}
  \item{selfemp}{factor. Is the individual self-employed?}
  \item{family}{family size.}
  \item{region}{factor indicating region.}
  \item{ethnicity}{factor indicating ethnicity: African-American, Caucasian, other.}
  \item{education}{factor indicating highest degree attained:
    no degree, GED (high school equivalent), high school, bachelor, master, PhD, other.}
}
}

\details{
This is a subset of the data used in Perry and Rosen (2004).
}

\source{
Online complements to Stock and Watson (2007). 

\url{http://wps.aw.com/aw_stock_ie_2/0,12040,3332253-,00.html}
}

\references{
Perry, C. and Rosen, H.S. (2004). \dQuote{The Self-Employed are Less Likely than
Wage-Earners to Have Health Insurance. So What?} in
Holtz-Eakin, D. and Rosen, H.S. (eds.), \emph{Entrepeneurship and
Public Policy}, MIT Press. 

Stock, J.H. and Watson, M.W. (2007). \emph{Introduction to Econometrics}, 2nd ed. Boston: Addison Wesley.
}

\seealso{\code{\link{StockWatson2007}}}

\examples{
data("HealthInsurance")
summary(HealthInsurance)
prop.table(xtabs(~ selfemp + insurance, data = HealthInsurance), 1)
}

\keyword{datasets}
