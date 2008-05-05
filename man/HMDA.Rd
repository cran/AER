\name{HMDA} 
\alias{HMDA} 
\title{Home Mortgage Disclosure Act Data} 

\description{Cross-section data on the Home Mortgage Disclosure Act (HMDA).}

\usage{data("HMDA")}

\format{
A data frame containing 2,380 observations on 14 variables.
\describe{
  \item{deny}{Factor. Was the mortgage denied?}
  \item{pirat}{Payments to income ratio.}    
  \item{hirat}{Housing expense to income ratio.}
  \item{lvrat}{Loan to value ratio.}
  \item{chist}{Factor. Credit history: consumer payments.}      
  \item{mhist}{Factor. Credit history: mortgage payments.}          
  \item{phist}{Factor. Public bad credit record?}    
  \item{unemp}{1989 Massachusetts unemployment rate in applicant's industry.}    
  \item{selfemp}{Factor. Is the individual self-employed?}  
  \item{insurance}{Factor. Was the individual denied mortgage insurance?}  
  \item{condomin}{Factor. Is the unit a condominium?}    
  \item{afam}{Factor. Is the individual African-American?}  
  \item{single}{Factor. Is the individual single?}    
  \item{hschool}{Factor. Does the individual have a high-school diploma?}    
}
}

\details{Only includes variables used by Stock and Watson (2007), 
some of which had to be generated from the raw data.
}

\source{
Online complements to Stock and Watson (2007). 

\url{http://wps.aw.com/aw_stock_ie_2}
}

\references{
Munnell, A. H., Tootell, G. M. B., Browne, L. E. and McEneaney, J. (1996). 
Mortgage Lending in Boston: Interpreting HMDA Data. \emph{American Economic Review}, \bold{86}, 25--53.

Stock, J. H. and Watson, M. W. (2007). \emph{Introduction to Econometrics}, 2nd ed. Boston: Addison Wesley.
}

\seealso{\code{\link{StockWatson2007}}}

\examples{
data("HMDA")

## Stock and Watson (2007)
## Equations 11.1, 11.3, 11.7, 11.8 and 11.10, pp. 387--395
fm1 <- lm(I(as.numeric(deny) - 1) ~ pirat, data = HMDA)
fm2 <- lm(I(as.numeric(deny) - 1) ~ pirat + afam, data = HMDA)
fm3 <- glm(deny ~ pirat, family = binomial(link = "probit"), data = HMDA)
fm4 <- glm(deny ~ pirat + afam, family = binomial(link = "probit"), data = HMDA)
fm5 <- glm(deny ~ pirat + afam, family = binomial(link = "logit"), data = HMDA)

## More examples can be found in:
## help("StockWatson2007")
}

\keyword{datasets}
