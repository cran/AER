\name{Affairs}
\alias{Affairs}
\title{Fair's Extramarital Affairs Data}

\description{
Infidelity data, known as Fair's Affairs. Cross-section data from a survey
conducted by Psychology Today in 1969.
}

\usage{data("Affairs")}

\format{
A data frame containing 601 observations on 9 variables.
\describe{
  \item{affairs}{numeric. How often engaged in extramarital sexual intercourse
    during the past year?}
  \item{gender}{factor indicating gender.}
  \item{age}{numeric variable coding age in years: \code{17.5} = under 20, \code{22} = 20--24,
    \code{27} = 25--29, \code{32} = 30--34, \code{37} = 35--39, \code{42} = 40--44,
    \code{47} = 45--49, \code{52} = 50--54, \code{57} = 55 or over.}
  \item{yearsmarried}{numeric variable coding number of years married: \code{0.125} = 3 months or less,
    \code{0.417} = 4--6 months, \code{0.75} = 6 months--1 year, \code{1.5} = 1--2 years,
    \code{4} = 3--5 years, \code{7} = 6--8 years, \code{10} = 9--11 years, \code{15} = 12 or more years.}
  \item{children}{factor. Are there children in the marriage?}
  \item{religiousness}{numeric variable coding religousness: \code{1} = anti, \code{2} = not at all,
    \code{3} = slightly, \code{4} = somewhat, \code{5} = very.}
  \item{education}{numeric variable coding level of education: \code{9} = grade school,
    \code{12} = high school graduate, \code{14} = some college, \code{16} = college graduate,
    \code{17} = some graduate work, \code{18} = master's degree, \code{20} = Ph.D., M.D., or
    other advanced degree.}
  \item{occupation}{numeric variable coding occupation according to Hollingshead classification
    (reverse numbering).}
  \item{rating}{numeric variable coding self rating of marriage: \code{1} = very unhappy,
    \code{2} = somewhat unhappy, \code{3} = average, \code{4} = happier than average, 
    \code{5} = very happy.}
}
}

\source{
Online complements to Greene (2003). Table F22.2.

\url{http://pages.stern.nyu.edu/~wgreene/Text/tables/tablelist5.htm}
}

\references{
Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Upper Saddle River, NJ: Prentice Hall.

Fair, R.C. (1978). A Theory of Extramarital Affairs. \emph{Journal of Political Economy},
\bold{86}, 45--61.
}

\seealso{\code{\link{Greene2003}}}

\examples{
data("Affairs")

## Greene (2003)
## Tab. 22.3 and 22.4
fm_ols <- lm(affairs ~ age + yearsmarried + religiousness + occupation + rating, data = Affairs)
fm_probit <- glm(I(affairs > 0) ~ age + yearsmarried + religiousness + occupation + rating,
  data = Affairs, family = binomial(link = "probit"))

fm_tobit <- tobit(affairs ~ age + yearsmarried + religiousness + occupation + rating, data = Affairs)
fm_tobit2 <- tobit(affairs ~ age + yearsmarried + religiousness + occupation + rating,
  right = 4, data = Affairs)

fm_pois <- glm(affairs ~ age + yearsmarried + religiousness + occupation + rating,
  data = Affairs, family = poisson)

library("MASS")
fm_nb <- glm.nb(affairs ~ age + yearsmarried + religiousness + occupation + rating, data = Affairs)

## Tab. 22.6
library("pscl")
fm_zip <- zeroinfl(affairs ~ age + yearsmarried + religiousness + occupation + rating | age + 
  yearsmarried + religiousness + occupation + rating, data = Affairs)
}

\keyword{datasets}
