\name{MarkPound}
\alias{MarkPound}

\title{DEM/GBP Exchange Rate Returns}

\description{
A daily time series of percentage returns of Deutsche mark/British pound (DEM/GBP)
exchange rates from 1984-01-03 through 1991-12-31.
}

\usage{data("MarkPound")}
\format{
A univariate time series of 1974 returns (exact dates unknown) for the DEM/GBP exchange rate.
}

\details{
Greene (2003, Table F11.1) rounded the series to six digits while eight digits are given in
Bollerslev and Ghysels (1996). Here, we provide the original data. Using \code{\link{round}}
a series can be produced that is virtually identical to that of Greene (2003) (except for
eight observations where a slightly different rounding arithmetic was used).
}

\source{
Journal of Business \& Economic Statistics Data Archive. 

\url{http://www.amstat.org/publications/jbes/upload/index.cfm?fuseaction=ViewArticles&pub=JBES&issue=96-2-APR}
}

\references{
Bollerslev, T., and Ghysels, E. (1996). Periodic Autoregressive Conditional Heteroskedasticity.  
\emph{Journal of Business \& Economic Statistics},
\bold{14}, 139--151.

Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Upper Saddle River, NJ: Prentice Hall.
}

\seealso{\code{\link{Greene2003}}, \code{\link{MarkDollar}}}

\examples{
## data as given by Greene (2003)
data("MarkPound")
mp <- round(MarkPound, digits = 6)

## Figure 11.3 in Greene (2003)
plot(mp)

## Example 11.8 in Greene (2003), Table 11.5
library("tseries")
mp_garch <- garch(mp, grad = "numerical")
summary(mp_garch)
logLik(mp_garch)  
## Greene (2003) also includes a constant and uses different
## standard errors (presumably computed from Hessian), here
## OPG standard errors are used. garchFit() in "fGarch"
## implements the approach used by Greene (2003).

## compare Errata to Greene (2003)
library("dynlm")
res <- residuals(dynlm(mp ~ 1))^2
mp_ols <- dynlm(res ~ L(res, 1:10))
summary(mp_ols)
logLik(mp_ols)
summary(mp_ols)$r.squared * length(residuals(mp_ols))
}

\keyword{datasets}
