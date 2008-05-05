\name{Medicaid1986}
\alias{Medicaid1986}
\title{Medicaid Utilization Data}
\description{
Cross-section data originating from the 1986 Medicaid Consumer Survey. The data comprise
two groups of Medicaid eligibles at two sites in California (Santa Barbara and Ventura 
counties): a group enrolled in a managed care demonstration program and a fee-for-service 
comparison group of non-enrollees.
}
\usage{data("Medicaid1986")}

\format{
A data frame containing 996 observations on 14 variables.
\describe{
  \item{visits}{Number of doctor visits.}
  \item{exposure}{Length of observation period for ambulatory care (days).}
  \item{children}{Total number of children in the household.}
  \item{age}{Age of the respondent.}
  \item{income}{Annual household income (average of income range in million USD).}
  \item{health1}{The first principal component (divided by 1000) of three
    health-status variables: functional limitations, acute conditions, and chronic
    conditions.}
  \item{health2}{The second principal component (divided by 1000) of three
    health-status variables: functional limitations, acute conditions, and chronic
    conditions.}
  \item{access}{Availability of health services (0 = low access, 1 = high access).}
  \item{married}{Factor. Is the individual married?}   
  \item{gender}{Factor indicating gender.}  
  \item{ethnicity}{Factor indicating ethnicity (\code{"cauc"} or \code{"other"}).}
  \item{school}{Number of years completed in school.}
  \item{enroll}{Factor. Is the individual enrolled in a demonstration program?}
  \item{program}{Factor indicating the managed care demonstration program:
    Aid to Families with Dependent Children (\code{"afdc"}) or
    non-institutionalized Supplementary Security Income (\code{"ssi"}).}
}
}

\source{
Journal of Applied Econometrics Data Archive. 

\url{http://qed.econ.queensu.ca/jae/1997-v12.3/gurmu/}
}

\references{
Gurmu, S. (1997). Semi-Parametric Estimation of Hurdle Regression Models 
with an Application to Medicaid Utilization.  \emph{Journal of Applied Econometrics},
\bold{12}, 225--242.
}

\examples{
## data and packages
data("Medicaid1986")
library("pscl")

## scale regressors
Medicaid1986$age2 <- Medicaid1986$age^2 / 100
Medicaid1986$school <- Medicaid1986$school / 10
Medicaid1986$income <- Medicaid1986$income / 10

## subsets
afdc <- subset(Medicaid1986, program == "afdc")[, c(1, 3:4, 15, 5:9, 11:13)]
ssi <- subset(Medicaid1986, program == "ssi")[, c(1, 3:4, 15, 5:13)]

## Gurmu (1997):
## Table VI., Poisson and negbin models
afdc_pois <- glm(visits ~ ., data = afdc, family = poisson)
summary(afdc_pois)
coeftest(afdc_pois, vcov = sandwich)

afdc_nb <- glm.nb(visits ~ ., data = afdc)
ssi_pois <- glm(visits ~ ., data = ssi, family = poisson)
ssi_nb <- glm.nb(visits ~ ., data = ssi)

## Table VII., Hurdle models (without semi-parametric effects)
afdc_hurdle <- hurdle(visits ~ . | . - access, data = afdc, dist = "negbin")
ssi_hurdle <- hurdle(visits ~ . | . - access, data = ssi, dist = "negbin")

## Table VIII., Observed and expected frequencies
round(cbind(
  Observed = table(afdc$visits)[1:8],
  Poisson = sapply(0:7, function(x) sum(dpois(x, fitted(afdc_pois)))),
  Negbin = sapply(0:7, function(x) sum(dnbinom(x, mu = fitted(afdc_nb), size = afdc_nb$theta))),
  Hurdle = colSums(predict(afdc_hurdle, type = "prob")[,1:8])
  )/nrow(afdc), digits = 3) * 100
round(cbind(
  Observed = table(ssi$visits)[1:8],
  Poisson = sapply(0:7, function(x) sum(dpois(x, fitted(ssi_pois)))),
  Negbin = sapply(0:7, function(x) sum(dnbinom(x, mu = fitted(ssi_nb), size = ssi_nb$theta))),
  Hurdle = colSums(predict(ssi_hurdle, type = "prob")[,1:8])
  )/nrow(ssi), digits = 3) * 100
}

\keyword{datasets}
