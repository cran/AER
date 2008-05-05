\name{CPS1985}
\alias{CPS1985}
\title{Determinants of Wages Data (CPS 1985)}

\description{
Cross-section data originating from the May 1985 Current Population Survey
by the US Census Bureau (random sample drawn for Berndt 1991).
}

\usage{data("CPS1985")}

\format{
A data frame containing 534 observations on 11 variables.
\describe{
  \item{wage}{Wage (in dollars per hour).}
  \item{education}{Number of years of education.}
  \item{experience}{Number of years of potential work experience
    (\code{age - education - 6}).}
  \item{age}{Age in years.}
  \item{ethnicity}{Factor with levels \code{"cauc"}, \code{"hispanic"},
    \code{"other"}.}
  \item{region}{Factor. Does the individual live in the South?}
  \item{gender}{Factor indicating gender.}
  \item{occupation}{Factor with levels \code{"worker"} (tradesperson or assembly line worker),
    \code{"technical"} (technical or professional worker), \code{"services"} (service
    worker), \code{"office"} (office and clerical worker), \code{"sales"} (sales worker),
    \code{"management"} (management and administration).}
  \item{sector}{Factor with levels \code{"manufacturing"} (manufacturing or mining),
    \code{"construction"}, \code{"other"}.}
  \item{union}{Factor. Does the individual work on a union job?}
  \item{married}{Factor. Is the individual married?}
}
}

\source{
StatLib. 

\url{http://lib.stat.cmu.edu/datasets/CPS_85_Wages}
}

\references{
Berndt, E.R. (1991). \emph{The Practice of Econometrics}. New York: Addison-Wesley.
}

\seealso{\code{\link{CPS1988}}, \code{\link{CPSSW}}}

\examples{
data("CPS1985")

## Berndt (1991)
## Exercise 2, p. 196
cps_2b <- lm(log(wage) ~ union + education, data = CPS1985)
cps_2c <- lm(log(wage) ~ -1 + union + education, data = CPS1985)

## Exercise 3, p. 198/199
cps_3a <- lm(log(wage) ~ education + experience + I(experience^2), data = CPS1985)
cps_3b <- lm(log(wage) ~ gender + education + experience + I(experience^2), data = CPS1985)
cps_3c <- lm(log(wage) ~ gender + married + education + experience + I(experience^2), data = CPS1985)
cps_3e <- lm(log(wage) ~ gender*married + education + experience + I(experience^2), data = CPS1985)

## Exercise 4, p. 199/200
cps_4a <- lm(log(wage) ~ gender + union + ethnicity + education + experience + I(experience^2), data = CPS1985)
cps_4c <- lm(log(wage) ~ gender + union + ethnicity + education * experience + I(experience^2), data = CPS1985)

## Exercise 6, p. 203
cps_6a <- lm(log(wage) ~ gender + union + ethnicity + education + experience + I(experience^2), data = CPS1985)
cps_6a_noeth <- lm(log(wage) ~ gender + union + education + experience + I(experience^2), data = CPS1985)
anova(cps_6a_noeth, cps_6a)

## Exercise 8, p. 208
cps_8a <- lm(log(wage) ~ gender + union + ethnicity + education + experience + I(experience^2), data = CPS1985)
summary(cps_8a)
coeftest(cps_8a, vcov = vcovHC(cps_8a, type = "HC0"))
}

\keyword{datasets}
