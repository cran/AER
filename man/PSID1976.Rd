\name{PSID1976}
\alias{PSID1976}

\title{Labor Force Participation Data}

\description{
Cross-section data originating from the 1976 Panel Study of Income Dynamics (PSID),
based on data for the previous year, 1975. 
}

\usage{data("PSID1976")}

\format{
A data frame containing 753 observations on 21 variables.
\describe{
  \item{participation}{Factor. Did the individual participate in the labor force in 1975?}
  \item{hours}{Wife's hours of work in 1975.}
  \item{youngkids}{Number of children less than 6 years old in household.}
  \item{oldkids}{Number of children between ages 6 and 18 in household.}
  \item{age}{Wife's age in years.}
  \item{education}{Wife's education in years.}
  \item{earnings}{Wife's average hourly earnings, in 1975 dollars.}
  \item{wage}{Wife's wage reported at the time of the 1976 interview 
    (not the same as the 1975 estimated wage). To use the subsample with this wage,
    one needs to select 1975 workers with \code{participation == "yes"}, then select only those
    women with non-zero wage. Only 325 women work in 1975 and have a non-zero wage in 1976.}
  \item{hhours}{Husband's hours worked in 1975.}
  \item{hage}{Husband's age in years.}
  \item{heducation}{Husband's education in years.}
  \item{hwage}{Husband's wage, in 1975 dollars.}
  \item{fincome}{Family income, in 1975 dollars. (This variable is used to construct the
    property income variable.)}
  \item{tax}{Marginal tax rate facing the wife, and is taken from published federal 
    tax tables (state and local income taxes are excluded). The taxable income on 
    which this tax rate is calculated includes Social Security, if applicable to wife.}
  \item{meducation}{Wife's mother's educational attainment, in years.}
  \item{feducation}{Wife's father's educational attainment, in years.}
  \item{unemp}{Unemployment rate in county of residence, in percentage points.
    (This taken from bracketed ranges.)}
  \item{city}{Factor. Does the individual live in a large city?}
  \item{experience}{Actual years of wife's previous labor market experience.}
  \item{college}{Factor. Did the individual attend college?}
  \item{hcollege}{Factor. Did the individual's husband attend college?}
}
}

\details{
This data set is also known under the name Mroz (1987) data.
Note that the variable \code{earnings} is called \code{wage} in many versions of this data set.
}

\source{
Online complements to Greene (2003). Table F4.1.

\url{http://pages.stern.nyu.edu/~wgreene/Text/tables/tablelist5.htm}
}

\references{
Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Upper Saddle River, NJ: Prentice Hall.

McCullough, B.D. (2004). Some Details of Nonlinear Estimation. In:
Altman, M., Gill, J., and McDonald, M.P.:  
\emph{Numerical Issues in Statistical Computing for the Social Scientist}. 
Hoboken, NJ: John Wiley, Ch. 8, 199--218.

Mroz, T.A. (1987). The Sensitivity of an Empirical Model of Married Women's Hours of Work to
Economic and Statistical Assumptions. \emph{Econometrica}, \bold{55}, 765--799.

Wooldridge, J.M. (2002). \emph{Econometric Analysis of Cross-Section and Panel Data}. 
Cambridge, MA: MIT Press.
}

\seealso{\code{\link{Greene2003}}, \code{\link{PSID1982}}}

\examples{
## data and transformations
data("PSID1976")
PSID1976$kids <- with(PSID1976, factor((youngkids + oldkids) > 0,
  levels = c(FALSE, TRUE), labels = c("no", "yes")))
PSID1976$nwincome <- with(PSID1976, (fincome - hours * earnings)/1000)
PSID1976$partnum <- as.numeric(PSID1976$participation) - 1

###################
## Greene (2003) ##
###################

## Example 4.1, Table 4.2
## (reproduced in Example 7.1, Table 7.1)
gr_lm <- lm(log(hours * earnings) ~ age + I(age^2) + education + kids,
  data = PSID1976, subset = participation == "yes")
summary(gr_lm)
vcov(gr_lm)

## Example 4.5
summary(gr_lm)
## or equivalently
gr_lm1 <- lm(log(hours * earnings) ~ 1, data = PSID1976, subset = participation == "yes")
anova(gr_lm1, gr_lm)

## Example 21.4, p. 681
gr_probit1 <- glm(participation ~ age + I(age^2) + fincome + education + kids,
  data = PSID1976, family = binomial(link = "probit") )
gr_probit2 <- glm(participation ~ age + I(age^2) + fincome + education,
  data = PSID1976, family = binomial(link = "probit"))
gr_probit3 <- glm(participation ~ kids/(age + I(age^2) + fincome + education),
  data = PSID1976, family = binomial(link = "probit"))
## LR test of all coefficients
lrtest(gr_probit1)
## Chow-type test
lrtest(gr_probit2, gr_probit3)
## equivalently:
anova(gr_probit2, gr_probit3, test = "Chisq")
## Table 21.3
summary(gr_probit1)

## Example 22.8, Table 22.7, p. 786
library("sampleSelection")
gr_2step <- selection(participation ~ age + I(age^2) + fincome + education + kids, 
  earnings ~ experience + I(experience^2) + education + city,
  data = PSID1976, method = "2step")
gr_ml <- selection(participation ~ age + I(age^2) + fincome + education + kids, 
  earnings ~ experience + I(experience^2) + education + city,
  data = PSID1976, method = "ml")
gr_ols <- lm(wage ~ experience + I(experience^2) + education + city,
  data = PSID1976, subset = participation == "yes")
## NOTE: huge discrepancies for ML and OLS, source unclear


#######################
## Wooldridge (2002) ##
#######################

## Table 15.1, p. 468
wl_lpm <- lm(partnum ~ nwincome + education + experience + I(experience^2) +
  age + youngkids + oldkids, data = PSID1976)
wl_logit <- glm(participation ~ nwincome + education + experience + I(experience^2) +
  age + youngkids + oldkids, family = binomial, data = PSID1976)
wl_probit <- glm(participation ~ nwincome + education + experience + I(experience^2) +
  age + youngkids + oldkids, family = binomial(link = "probit"), data = PSID1976)
## (same as Altman et al.)

## convenience functions
pseudoR2 <- function(obj) 1 - as.vector(logLik(obj)/logLik(update(obj, . ~ 1)))
misclass <- function(obj) 1 - sum(diag(prop.table(table(
  model.response(model.frame(obj)), round(fitted(obj))))))

coeftest(wl_logit)
logLik(wl_logit)
misclass(wl_logit)
pseudoR2(wl_logit)

coeftest(wl_probit)
logLik(wl_probit)
misclass(wl_probit)
pseudoR2(wl_probit)

## Table 16.2, p. 528
form <- hours ~ nwincome + education + experience + I(experience^2) + age + youngkids + oldkids 
wl_ols <- lm(form, data = PSID1976)
wl_tobit <- tobit(form, data = PSID1976)
summary(wl_ols)
summary(wl_tobit)


#######################
## McCullough (2004) ##
#######################

## p. 203
mc_probit <- glm(participation ~ nwincome + education + experience + I(experience^2) +
  age + youngkids + oldkids, family = binomial(link = "probit"), data = PSID1976)
mc_tobit <- tobit(hours ~ nwincome + education + experience + I(experience^2) + age +
  youngkids + oldkids, data = PSID1976)
coeftest(mc_probit)
coeftest(mc_tobit)
coeftest(mc_tobit, vcov = vcovOPG)
}

\keyword{datasets}
