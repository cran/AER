\name{NMES1988}
\alias{NMES1988}
\title{Demand for Medical Care in NMES 1988}
\description{
Cross-section data originating from the US National Medical Expenditure Survey (NMES) 
conducted in 1987 and 1988. The NMES is based upon a representative, national 
probability sample of the civilian non-institutionalized population and individuals 
admitted to long-term care facilities during 1987. The data are a subsample of 
individuals ages 66 and over all of whom are covered by Medicare 
(a public insurance program providing substantial protection against health-care costs).
}
\usage{data("NMES1988")}

\format{
A data frame containing 4,406 observations on 19 variables.
\describe{
  \item{visits}{Number of physician office visits.}
  \item{nvisits}{Number of non-physician office visits.}
  \item{ovisits}{Number of physician hospital outpatient visits.}
  \item{novisits}{Number of non-physician hospital outpatient visits.}
  \item{emergency}{Emergency room visits.}
  \item{hospital}{Number of hospital stays.}
  \item{health}{Factor indicating self-perceived health status, levels are
    \code{"poor"}, \code{"average"} (reference category), \code{"excellent"}.}
  \item{chronic}{Number of chronic conditions.}
  \item{adl}{Factor indicating whether the individual has a condition that
    limits activities of daily living (\code{"limited"}) or not (\code{"normal"}).}
  \item{region}{Factor indicating region, levels are \code{northeast},
    \code{midwest}, \code{west}, \code{other} (reference category).}  
  \item{age}{Age in years (divided by 10).}
  \item{afam}{Factor. Is the individual African-American?}
  \item{gender}{Factor indicating gender.}
  \item{married}{Factor. is the individual married?}
  \item{school}{Number of years of education.}
  \item{income}{Family income in USD 10,000.}
  \item{employed}{Factor. Is the individual employed?}
  \item{insurance}{Factor. Is the individual covered by private insurance?}
  \item{medicaid}{Factor. Is the individual covered by Medicaid?}
}
}

\source{
Journal of Applied Econometrics Data Archive for Deb and Trivedi (1997). 

\url{http://www.econ.queensu.ca/jae/1997-v12.3/deb-trivedi/}
}

\references{
Cameron, A.C. and Trivedi, P.K. (1998). \emph{Regression Analysis of Count Data}. 
Cambridge: Cambridge University Press.

Deb, P., and Trivedi, P.K. (1997). Demand for Medical Care by the
Elderly: A Finite Mixture Approach. \emph{Journal of Applied Econometrics},
\bold{12}, 313--336.

Zeileis, A., Kleiber, C., and Jackman, S. (2008). Regression Models
for Count Data in R. \emph{Journal of Statistical Software},
forthcoming.
}

\seealso{\code{\link{CameronTrivedi1998}}}

\examples{
## package
library("pscl")

## select variables for analysis
data("NMES1988")
nmes <- NMES1988[, c(1, 6:8, 13, 15, 18)]

## dependent variable
hist(nmes$visits, breaks = 0:(max(nmes$visits)+1) - 0.5)
plot(table(nmes$visits))

## convenience transformations for exploratory graphics
clog <- function(x) log(x + 0.5)
cfac <- function(x, breaks = NULL) {
  if(is.null(breaks)) breaks <- unique(quantile(x, 0:10/10))
  x <- cut(x, breaks, include.lowest = TRUE, right = FALSE)
  levels(x) <- paste(breaks[-length(breaks)], ifelse(diff(breaks) > 1,
    c(paste("-", breaks[-c(1, length(breaks))] - 1, sep = ""), "+"), ""), sep = "")
  return(x)
}

## bivariate visualization
par(mfrow = c(3, 2))
plot(clog(visits) ~ health, data = nmes, varwidth = TRUE)
plot(clog(visits) ~ cfac(chronic), data = nmes)
plot(clog(visits) ~ insurance, data = nmes, varwidth = TRUE)
plot(clog(visits) ~ cfac(hospital, c(0:2, 8)), data = nmes)
plot(clog(visits) ~ gender, data = nmes, varwidth = TRUE)
plot(cfac(visits, c(0:2, 4, 6, 10, 100)) ~ school, data = nmes, breaks = 9)
par(mfrow = c(1, 1))

## Poisson regression
nmes_pois <- glm(visits ~ ., data = nmes, family = poisson)
summary(nmes_pois)

## LM test for overdispersion
dispersiontest(nmes_pois)
dispersiontest(nmes_pois, trafo = 2)

## sandwich covariance matrix
coeftest(nmes_pois, vcov = sandwich)

## quasipoisson model
nmes_qpois <- glm(visits ~ ., data = nmes, family = quasipoisson)

## NegBin regression
nmes_nb <- glm.nb(visits ~ ., data = nmes)

## hurdle regression
nmes_hurdle <- hurdle(visits ~ . | hospital + chronic + insurance + school + gender,
  data = nmes, dist = "negbin")

## zero-inflated regression model
nmes_zinb <- zeroinfl(visits ~ . | hospital + chronic + insurance + school + gender,
  data = nmes, dist = "negbin")

## compare estimated coefficients
fm <- list("ML-Pois" = nmes_pois, "Quasi-Pois" = nmes_qpois, "NB" = nmes_nb,
  "Hurdle-NB" = nmes_hurdle, "ZINB" = nmes_zinb)
round(sapply(fm, function(x) coef(x)[1:8]), digits = 3)

## associated standard errors
round(cbind("ML-Pois" = sqrt(diag(vcov(nmes_pois))),
  "Adj-Pois" = sqrt(diag(sandwich(nmes_pois))),
  sapply(fm[-1], function(x) sqrt(diag(vcov(x)))[1:8])),
  digits = 3)

## log-likelihoods and number of estimated parameters
rbind(logLik = sapply(fm, function(x) round(logLik(x), digits = 0)),
  Df = sapply(fm, function(x) attr(logLik(x), "df")))

## predicted number of zeros
round(c("Obs" = sum(dt$visits < 1),
  "ML-Pois" = sum(dpois(0, fitted(nmes_pois))),
  "Adj-Pois" = NA,
  "Quasi-Pois" = NA,
  "NB" = sum(dnbinom(0, mu = fitted(nmes_nb), size = nmes_nb$theta)),
  "NB-Hurdle" = sum(predict(nmes_hurdle, type = "prob")[,1]),
  "ZINB" = sum(predict(nmes_zinb, type = "prob")[,1])))

## coefficients of zero-augmentation models
t(sapply(fm[4:5], function(x) round(x$coefficients$zero, digits = 3)))
}

\keyword{datasets}
