\name{NaturalGas}
\alias{NaturalGas}

\title{Natural Gas Data}

\description{
Panel data originating from 6 US states over the period 1967--1989.
}

\usage{data("NaturalGas")}

\format{
A data frame containing 138 observations on 10 variables.
\describe{
  \item{state}{factor. State abbreviation.}
  \item{statecode}{factor. State Code.}
  \item{year}{factor coding year.}
  \item{consumption}{Consumption of natural gas by the residential sector.}
  \item{price}{Price of natural gas}
  \item{eprice}{Price of electricity.}
  \item{oprice}{Price of distillate fuel oil.}
  \item{lprice}{Price of liquefied petroleum gas.}
  \item{heating}{Heating degree days.}
  \item{income}{Real per-capita personal income.}
}
}

\source{
The data are from Baltagi (2002) and available at

\url{http://www.springeronline.com/sgw/cda/frontpage/0,10735,4-165-2-107420-0,00.html}
}

\references{
Baltagi, B.H. (2002). \emph{Econometrics}, 3rd ed. Berlin, Springer. 
}

\seealso{\code{\link{Baltagi2002}}}

\examples{
data("NaturalGas")
summary(NaturalGas)
}

\keyword{datasets}
