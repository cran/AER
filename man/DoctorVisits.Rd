\name{DoctorVisits}
\alias{DoctorVisits}

\title{Australian Health Service Utilization Data}

\description{
Cross-section data originating from the 1977--1978 Australian Health Survey. 
}

\usage{data("DoctorVisits")}

\format{
A data frame containing 5,190 observations on 12 variables.
\describe{
  \item{visits}{Number of doctor visits in past 2 weeks.}
  \item{gender}{Factor indicating gender.}
  \item{age}{Age in years divided by 100.}
  \item{income}{Annual income in tens of thousands of dollars.}
  \item{illness}{Number of illnesses in past 2 weeks.}
  \item{reduced}{Number of days of reduced activity in past 2 weeks due to illness or injury.}
  \item{health}{General health questionnaire score using Goldberg's method.}
  \item{private}{Factor. Does the individual have private health insurance?}
  \item{freepoor}{Factor. Does the individual have free government health insurance due to low income?}   
  \item{freerepat}{Factor. Does the individual have free government health insurance due to old age, 
          disability or veteran status?}  
  \item{nchronic}{Factor. Is there a chronic condition not limiting activity?}
  \item{lchronic}{Factor. Is there a chronic condition limiting activity?}
}
}

\source{
Journal of Applied Econometrics Data Archive. 

\url{http://www.econ.queensu.ca/jae/1997-v12.3/mullahy/}
}

\references{
Cameron, A.C. and Trivedi, P.K. (1986). Econometric Models Based on Count Data: 
Comparisons and Applications of Some Estimators and Tests.
\emph{Journal of Applied Econometrics},
\bold{1}, 29--53.

Cameron, A.C. and Trivedi, P.K. (1998). \emph{Regression Analysis of Count Data}. 
Cambridge: Cambridge University Press.

Mullahy, J. (1997). Heterogeneity, Excess Zeros, and the Structure of Count Data Models.  
\emph{Journal of Applied Econometrics},
\bold{12}, 337--350.
}

\seealso{\code{\link{CameronTrivedi1998}}}

\examples{
data("DoctorVisits")
library("MASS")

## Cameron and Trivedi (1986), Table III, col. (1)
dv_lm <- lm(visits ~ . + I(age^2), data = DoctorVisits)
summary(dv_lm)

## Cameron and Trivedi (1998), Table 3.3 
dv_pois <- glm(visits ~ . + I(age^2), data = DoctorVisits, family = poisson)
summary(dv_pois)                  ## MLH standard errors
coeftest(dv_pois, vcov = vcovOPG) ## MLOP standard errors
logLik(dv_pois)
## standard errors denoted RS ("unspecified omega robust sandwich estimate")
coeftest(dv_pois, vcov = sandwich)

## Cameron and Trivedi (1986), Table III, col. (4)
dv_nb <- glm.nb(visits ~ . + I(age^2), data = DoctorVisits)
summary(dv_nb)
logLik(dv_nb)
}

\keyword{datasets}
