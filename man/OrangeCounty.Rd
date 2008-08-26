\name{OrangeCounty}
\alias{OrangeCounty}

\title{Orange County Employment}

\description{
Quarterly time series data on employment in Orange county, 1965--1983.
}

\usage{data("OrangeCounty")}

\format{
A quarterly multiple time series from 1965 to 1983 with 2 variables.
\describe{
  \item{employment}{Quarterly employment in Orange county.}
  \item{gnp}{Quarterly real GNP.}
}
}

\source{
The data is from Baltagi (2002) and available at

\url{http://www.springeronline.com/sgw/cda/frontpage/0,10735,4-165-2-107420-0,00.html}
}

\references{
Baltagi, B.H. (2002). \emph{Econometrics}, 3rd ed. Berlin, Springer. 

}

\seealso{\code{\link{Baltagi2002}}}

\examples{
data("OrangeCounty")
plot(OrangeCounty)
}

\keyword{datasets}
