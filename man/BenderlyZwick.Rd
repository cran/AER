\name{BenderlyZwick}
\alias{BenderlyZwick}

\title{Benderly and Zwick Data: Inflation, Growth and Stock Returns}

\description{
Time-series data, 1952--1982.  
}
\usage{data("BenderlyZwick")}
\format{
An annual multiple time series from 1952 to 1982 with 3 variables.
\describe{
  \item{returns}{real annual returns on stocks, measured using
    the Ibbotson-Sinquefeld data base.}
  \item{growth}{annual growth rate of output, measured by real GNP
    (from the given year to the next year).}
  \item{inflation}{inflation rate, measured as growth of price
    rate (from December of the previous year to December of the present year).}
  \item{growth2}{annual growth rate of real GNP as given by Baltagi.}
  \item{inflation2}{inflation rate as given by Baltagi}
}
}

\source{
The first three columns of the data is from Table 1 in Benderly and Zwick (1985).
The remaining columns are taken from the online complements of Baltagi (2002).
The first column is identical in both sources, the other two variables differ
in their numeric values and additionally the growth series seems to be lagged
differently. Baltagi (2002) states Lott and Ray (1992) as the source for his
version of the data set which is available from

\url{http://www.springeronline.com/sgw/cda/frontpage/0,10735,4-165-2-107420-0,00.html}
}

\references{
Baltagi, B.H. (2002). \emph{Econometrics}, 3rd ed. Berlin, Springer. 

Benderly, J., and Zwick, B. (1985). Inflation, Real Balances, Output and Real Stock
Returns. \emph{American Economic Review}, \bold{75}, 1115--1123.

Lott, W.F., and Ray, S.C. (1992). \emph{Applied Econometrics: Problems with Data Sets}.
New York: The Dryden Press.

Zaman, A., Rousseeuw, P.J., and  Orhan, M. (2001). Econometric Applications of
High-Breakdown Robust Regression Rechniques. \emph{Economics Letters}, \bold{71}, 1--8.
}

\seealso{\code{\link{Baltagi2002}}}

\examples{
data("BenderlyZwick")
plot(BenderlyZwick)

## Benderly and Zwick (1985), p. 1116
library("dynlm")
bz_ols <- dynlm(returns ~ growth + inflation,
  data = BenderlyZwick/100, start = 1956, end = 1981)
summary(bz_ols)

## Zaman, Rousseeuw and Orhan (2001)
## use larger period, without scaling
bz_ols2 <- dynlm(returns ~ growth + inflation,
  data = BenderlyZwick, start = 1954, end = 1981)
summary(bz_ols2)
}

\keyword{datasets}
