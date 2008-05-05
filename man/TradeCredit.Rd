\name{TradeCredit}
\alias{TradeCredit}
\title{Trade Credit Data}

\description{
Macroeconomic time-series data from 1946 to 1966 on trade credit and the money market.
}

\usage{data("TradeCredit")}

\format{
An annual multiple time series from 1946 to 1966 on 7 variables.
\describe{
\item{trade}{Nominal total trade money.}
\item{reserve}{Nominal effective reserve money.}
\item{gnp}{GNP in current dollars.}
\item{utilization}{Degree of market utilization.}
\item{interest}{Short-term rate of interest.}
\item{size}{Mean real size of the representative economic unit (1939 = 100).}
\item{price}{GNP price deflator (1958 = 100).}
}
}

\source{
The data is from Baltagi (2002) and available at

\url{http://www.springeronline.com/sgw/cda/frontpage/0,10735,4-165-2-107420-0,00.html}
}

\references{
Baltagi, B.H. (2002). \emph{Econometrics}, 3rd ed. Berlin, Springer. 

Laffer, A.B. (1970). Trade Credit and the Money Market. \emph{Journal of Political Economy},
\bold{78}, 239--267.

}

\seealso{\code{\link{Baltagi2002}}}

\examples{
data("TradeCredit")
plot(TradeCredit)
}

\keyword{datasets}
