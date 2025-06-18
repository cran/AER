\name{PSID7682}
\alias{PSID7682}

\title{PSID Earnings Panel Data (1976--1982)}

\description{
Panel data on earnings of 595 individuals for the years 1976--1982,
originating from the Panel Study of Income Dynamics.  
}

\usage{data("PSID7682")}

\format{
A data frame containing 7 annual observations on 12 variables for 595 individuals. 
\describe{
  \item{experience}{Years of full-time work experience.}
  \item{weeks}{Weeks worked.}
  \item{occupation}{factor. Is the individual a white-collar (\code{"white"})
    or blue-collar (\code{"blue"}) worker?}
  \item{industry}{factor. Does the individual work in a manufacturing industry?}
  \item{south}{factor. Does the individual reside in the South?}
  \item{smsa}{factor. Does the individual reside in a SMSA (standard
    metropolitan statistical area)?}
  \item{married}{factor. Is the individual married?}
  \item{gender}{factor indicating gender.}
  \item{union}{factor. Is the individual's wage set by a union contract?}
  \item{education}{Years of education.}
  \item{ethnicity}{factor indicating ethnicity.
    Is the individual African-American (\code{"afam"}) or not (\code{"other"})?}
  \item{wage}{Wage.}
  \item{year}{factor indicating year.}
  \item{id}{factor indicating individual subject ID.}
}
}

\details{
The data were originally analyzed by Cornwell and Rupert (1988) and employed
for assessing various instrumental-variable estimators for panel models
(including the Hausman-Taylor model). Baltagi and Khanti-Akom (1990) reanalyzed
the data, made corrections to the data and also suggest modeling with a different
set of instruments.

\code{PSID7682} is the version of the data as provided by Baltagi (2005),
or Greene (2008).

Baltagi (2002) just uses the cross-section for the year 1982,
i.e., \code{subset(PSID7682, year == "1982")}. This is also available as
a standalone data set \code{\link{PSID1982}} because it was included
in \pkg{AER} prior to the availability of the full \code{PSID7682} panel
version.
}

\source{
Online complements to Baltagi (2005).

\url{https://www.wiley.com/legacy/wileychi/baltagi3e/data_sets.html}

Also provided in the online complements to Greene (2008), Table F9.1.

\url{https://pages.stern.nyu.edu/~wgreene/Text/Edition6/tablelist6.htm}
}

\references{
Baltagi, B.H., and Khanti-Akom, S. (1990). On Efficient Estimation with
Panel Data: An Empirical Comparison of Instrumental Variables Estimators.
\emph{Journal of Applied Econometrics}, \bold{5}, 401--406.

Baltagi, B.H. (2001). \emph{Econometric Analysis of Panel Data}, 2nd ed. Chichester, UK: John Wiley.

Baltagi, B.H. (2002). \emph{Econometrics}, 3rd ed. Berlin, Springer. 

Baltagi, B.H. (2005). \emph{Econometric Analysis of Panel Data}, 3rd ed. Chichester, UK: John Wiley.

Cornwell, C., and Rupert, P. (1988). Efficient Estimation with Panel Data: 
An Empirical Comparison of Instrumental Variables Estimators. 
\emph{Journal of Applied Econometrics}, \bold{3}, 149--155.

Greene, W.H. (2008). \emph{Econometric Analysis}, 6th ed. Upper Saddle River, NJ: Prentice Hall.
}

\seealso{\code{\link{PSID1982}}, \code{\link{Baltagi2002}}}

\examples{
\dontshow{ if(!requireNamespace("plm")) {
  if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
    stop("not all packages required for the example are installed")
  } else q() }}
data("PSID7682")

library("plm")
psid <- pdata.frame(PSID7682, c("id", "year"))

## Baltagi & Khanti-Akom, Table I, column "HT"
## original Cornwell & Rupert choice of exogenous variables
psid_ht1 <- plm(log(wage) ~ weeks + south + smsa + married +
  experience + I(experience^2) + occupation + industry + union + gender + ethnicity + education |
  weeks + south + smsa + married + gender + ethnicity,
  data = psid, model = "ht")

## Baltagi & Khanti-Akom, Table II, column "HT"
## alternative choice of exogenous variables
psid_ht2 <- plm(log(wage) ~ occupation + south + smsa + industry +
  experience + I(experience^2) + weeks + married + union + gender + ethnicity + education |
  occupation + south + smsa + industry + gender + ethnicity,
  data = psid, model = "ht")

## Baltagi & Khanti-Akom, Table III, column "HT"
## original choice of exogenous variables + time dummies
## (see also Baltagi, 2001, Table 7.1)
psid$time <- psid$year
psid_ht3 <- plm(log(wage) ~ weeks + south + smsa + married + experience + I(experience^2) +
  occupation + industry + union + gender + ethnicity + education + time |
  weeks + south + smsa + married + gender + ethnicity + time,
  data = psid, model = "ht")
}

\keyword{datasets}
