\name{RecreationDemand}
\alias{RecreationDemand}
\title{Recreation Demand Data}
\description{
Cross-section data on the number of recreational boating trips to Lake Somerville, Texas, in 1980,  
based on a survey administered to 2,000 registered leisure boat owners in 23 counties in eastern Texas. 
}
\usage{data("RecreationDemand")}

\format{
A data frame containing 659 observations on 8 variables.
\describe{
  \item{trips}{Number of recreational boating trips.}
  \item{quality}{Facility's subjective quality ranking on a scale of 1 to 5.}
  \item{ski}{factor. Was the individual engaged in water-skiing at the lake?}
  \item{income}{Annual household income of the respondent (in 1,000 USD).}
  \item{userfee}{factor. Did the individual pay an annual user fee at Lake Somerville?}
  \item{costC}{Expenditure when visiting Lake Conroe (in USD).}
  \item{costS}{Expenditure when visiting Lake Somerville (in USD).}
  \item{costH}{Expenditure when visiting Lake Houston (in USD).}
}
}


\details{
According to the original source (Seller, Stoll and Chavas, 1985, p. 168), the quality rating is on a scale from 1 to 5 and gives 0 for those who had not visited the lake. This explains the remarkably low mean for this variable, but also suggests that its treatment in various more recent publications is far from ideal. For consistency with other sources we handle the variable as a numerical variable, including the zeros.   
}

\source{
Journal of Business \& Economic Statistics Data Archive. 

\url{http://www.amstat.org/publications/jbes/upload/index.cfm?fuseaction=ViewArticles&pub=JBES&issue=96-4-OCT}
}

\references{
Cameron, A.C. and Trivedi, P.K. (1998). \emph{Regression Analysis of Count Data}. 
Cambridge: Cambridge University Press.

Gurmu, S. and Trivedi, P.K. (1996). Excess Zeros in Count Models for Recreational Trips.  
\emph{Journal of Business \& Economic Statistics},
\bold{14}, 469--477.

Ozuna, T. and Gomez, I.A. (1995). Specification and Testing of Count Data Recreation Demand Functions.
\emph{Empirical Economics}, \bold{20}, 543--550.

Seller, C., Stoll, J.R. and Chavas, J.-P. (1985). Validation og Empirical Measures of Welfare Change: A Comparison of Nonmarket Techniques. \emph{Land Economics}, \bold{61}, 156--175.
}

\seealso{\code{\link{CameronTrivedi1998}}}

\examples{
data("RecreationDemand")

## Poisson model:
## Cameron and Trivedi (1998), Table 6.11
## Ozuna and Gomez (1995), Table 2, col. 3
fm_pois <- glm(trips ~ ., data = RecreationDemand, family = poisson)
summary(fm_pois)
logLik(fm_pois)
coeftest(fm_pois, vcov = sandwich)

## Negbin model:
## Cameron and Trivedi (1998), Table 6.11
## Ozuna and Gomez (1995), Table 2, col. 5
library("MASS")
fm_nb <- glm.nb(trips ~ ., data = RecreationDemand)
coeftest(fm_nb, vcov = vcovOPG)

## ZIP model:
## Cameron and Trivedi (1998), Table 6.11
library("pscl")
fm_zip <- zeroinfl(trips ~  . | quality + income, data = RecreationDemand)
summary(fm_zip)

## Hurdle models
## Cameron and Trivedi (1998), Table 6.13
## poisson-poisson
fm_hp <- hurdle(trips ~ ., data = RecreationDemand, dist = "poisson", zero = "poisson")
## negbin-negbin
fm_hnb <- hurdle(trips ~ ., data = RecreationDemand, dist = "negbin", zero = "negbin")
## binom-negbin == geo-negbin
fm_hgnb <- hurdle(trips ~ ., data = RecreationDemand, dist = "negbin")

## Note: quasi-complete separation
with(RecreationDemand, table(trips > 0, userfee))
}

\keyword{datasets}
