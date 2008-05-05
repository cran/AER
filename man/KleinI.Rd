\name{KleinI}
\alias{KleinI}

\title{Klein Model I}

\description{
Klein's Model I for the US economy.
}

\usage{data("KleinI")}

\format{
An annual multiple time series from 1920 to 1941 with 9 variables.
 \describe{
    \item{consumption}{Consumption.}
    \item{cprofits}{Corporate profits.}
    \item{pwage}{Private wage bill.}
    \item{invest}{Investment.}
    \item{capital}{Previous year's capital stock.}
    \item{gnp}{Gross national product.}    
    \item{gwage}{Government wage bill.}
    \item{gexpenditure}{Government spending.}
    \item{taxes}{Taxes.}
 }
}

\source{
Online complements to Greene (2003). Table F15.1.

\url{http://pages.stern.nyu.edu/~wgreene/Text/tables/tablelist5.htm}
}

\references{
Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Upper Saddle River, NJ: Prentice Hall.

Klein, L. (1950). \emph{Economic Fluctuations in the United States, 1921--1941}. New York: John Wiley.

Maddala, G.S. (1977). \emph{Econometrics}. New York: McGraw-Hill.
}

\seealso{\code{\link{Greene2003}}}

\examples{
data("KleinI", package = "AER")
plot(KleinI)

## Greene (2003), Tab. 15.3, OLS
library("dynlm")
fm_cons <- dynlm(consumption ~ cprofits + L(cprofits) + I(pwage + gwage), data = KleinI)
fm_inv <- dynlm(invest ~ cprofits + L(cprofits) + capital, data = KleinI)
fm_pwage <- dynlm(pwage ~ gnp + L(gnp) + I(time(gnp) - 1931), data = KleinI)
summary(fm_cons)
summary(fm_inv)
summary(fm_pwage)

## More examples can be found in:
## help("Greene2003")
}

\keyword{datasets}
