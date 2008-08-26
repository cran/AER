\name{USMacroB}
\alias{USMacroB}
\title{US Macroeconomic Data (1959--1995, Baltagi)}

\description{
Time series data on 3 US macroeconomic variables for 1959--1995,
extracted from the Citibank data base.
}

\usage{data("USMacroB")}

\format{
A quarterly multiple time series from 1959(1) to 1995(2) with 3 variables.
\describe{
  \item{gnp}{Gross national product.}
  \item{mbase}{Average of the seasonally adjusted monetary base.}
  \item{tbill}{Average of 3 month treasury-bill rate (per annum).}
}
}

\source{
The data is from Baltagi (2002) and available at

\url{http://www.springeronline.com/sgw/cda/frontpage/0,10735,4-165-2-107420-0,00.html}
}

\references{
Baltagi, B.H. (2002). \emph{Econometrics}, 3rd ed. Berlin, Springer. 

}

\seealso{\code{\link{Baltagi2002}}, \code{\link{USMacroSW}}, \code{\link{USMacroSWQ}},
  \code{\link{USMacroSWM}}, \code{\link{USMacroG}}}

\examples{
data("USMacroB")
plot(USMacroB)
}

\keyword{datasets}
