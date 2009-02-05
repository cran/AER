\name{GrowthDJ}
\alias{GrowthDJ}

\encoding{UTF-8}

\title{Determinants of Economic Growth}

\description{Growth regression data as provided by Durlauf & Johnson (1995).}

\usage{data("GrowthDJ")}

\format{
A data frame containing 121 observations on 10 variables.
\describe{
 \item{oil}{factor. Is the country an oil-producing country?}
 \item{inter}{factor. Does the country have better quality data?}
 \item{oecd}{factor. Is the country a member of the OECD?}
 \item{gdp60}{Per capita GDP in 1960.}
 \item{gdp85}{Per capita GDP in 1985.}
 \item{gdpgrowth}{Average growth rate of per capita GDP from 1960 to 1985 (in percent).}
 \item{popgrowth}{Average growth rate of working-age population 1960 to 1985 (in percent).}
 \item{invest}{Average ratio of investment (including Government Investment) to GDP 
   from 1960 to 1985 (in percent).}
 \item{school}{Average fraction of working-age population enrolled in secondary school
   from 1960 to 1985 (in percent).}   
 \item{literacy60}{Fraction of the population over 15 years old that is able to read 
   and write in 1960 (in percent).}  
}
}

\details{
The data are derived from the Penn World Table 4.0 and are given in Mankiw, Romer and Weil (1992),
except \code{literacy60} that is from the World Bank's World Development Report.
}

\source{
Journal of Applied Econometrics Data Archive. 

\url{http://www.econ.queensu.ca/jae/1995-v10.4/durlauf-johnson/}
}

\references{
Durlauf, S.N., and Johnson, P.A. (1995). Multiple Regimes and Cross-Country
Growth Behavior. \emph{Journal of Applied Econometrics}, \bold{10}, 365--384. 

Koenker, R., and Zeileis, A. (2007). \dQuote{Reproducible Econometric Research
(A Critical Review of the State of the Art).} Report 60, Department of Statistics
and Mathematics, Wirtschaftsuniversität Wien, Research Report Series.
\url{http://epub.wu-wien.ac.at/}.

Mankiw, N.G, Romer, D., and Weil, D.N. (1992). A Contribution to the Empirics of Economic Growth. 
\emph{Quarterly Journal of Economics}, \bold{107}, 407--437. 

Masanjala, W.H., and Papageorgiou, C. (2004). The Solow Model with CES
Technology: Nonlinearities and Parameter Heterogeneity. \emph{Journal of
Applied Econometrics}, \bold{19}, 171--201.
}

\seealso{\code{\link{OECDGrowth}}, \code{\link{GrowthSW}}}

\examples{
## data for non-oil-producing countries
data("GrowthDJ")
dj <- subset(GrowthDJ, oil == "no")
## Different scalings have been used by different authors,
## different types of standard errors, etc.,
## see Koenker & Zeileis (2007) for an overview

## Durlauf & Johnson (1995), Table II
mrw_model <- I(log(gdp85) - log(gdp60)) ~ log(gdp60) +
  log(invest/100) + log(popgrowth/100 + 0.05) + log(school/100)
dj_mrw <- lm(mrw_model, data = dj)
coeftest(dj_mrw) 

dj_model <- I(log(gdp85) - log(gdp60)) ~ log(gdp60) +
  log(invest) + log(popgrowth/100 + 0.05) + log(school)
dj_sub1 <- lm(dj_model, data = dj, subset = gdp60 < 1800 & literacy60 < 50)
coeftest(dj_sub1, vcov = sandwich)

dj_sub2 <- lm(dj_model, data = dj, subset = gdp60 >= 1800 & literacy60 >= 50)
coeftest(dj_sub2, vcov = sandwich)
}

\keyword{datasets}
