\name{MarkDollar}
\alias{MarkDollar}

\title{DEM/USD Exchange Rate Returns}

\description{
A time series of intra-day percentage returns of Deutsche mark/US dollar (DEM/USD)
exchange rates, consisting of two observations per day from 1992-10-01 through 1993-09-29.
}

\usage{data("MarkDollar")}

\format{
A univariate time series of 518 returns (exact dates unknown) for the DEM/USD exchange rate.
}

\source{
Journal of Business & Economic Statistics Data Archive. 

\verb{http://www.amstat.org/publications/jbes/upload/index.cfm?fuseaction=ViewArticles&pub=JBES&issue=96-2-APR}
}

\references{
Bollerslev, T., and Ghysels, E. (1996). Periodic Autoregressive Conditional Heteroskedasticity.  
\emph{Journal of Business & Economic Statistics},
\bold{14}, 139--151.
}

\seealso{\code{\link{MarkPound}}}

\examples{
\dontshow{ if(!requireNamespace("tseries")) {
  if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
    stop("not all packages required for the example are installed")
  } else q() }}
library("tseries")
data("MarkDollar")

## GARCH(1,1)
fm <- garch(MarkDollar, grad = "numerical")
summary(fm)
logLik(fm)  
}

\keyword{datasets}
