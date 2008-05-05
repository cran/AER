\name{CartelStability}
\alias{CartelStability}

\title{CartelStability}

\description{
Weekly observations on prices and other factors from 1880--1886, for a total of 326 weeks.
}

\usage{data("CartelStability")}

\format{
A data frame containing 328 observations on 5 variables.
\describe{
  \item{price}{weekly index of price of shipping a ton of grain by rail.}
  \item{cartel}{factor. Is a railroad cartel operative?}
  \item{quantity}{total tonnage of grain shipped in the week.}
  \item{season}{factor indicating season of year. To match the weekly data,
    the calendar has been divided into 13 periods, each approximately 4 weeks long.}
  \item{ice}{factor. Are the Great Lakes innavigable because of ice?}
}
}

\source{
Online complements to Stock and Watson (2007). 

\url{http://wps.aw.com/aw_stock_ie_2/0,12040,3332253-,00.html}
}

\references{
Porter, R. H. (1983). A Study of Cartel Stability: The Joint Executive Committee, 1880--1886. \emph{The Bell Journal of Economics}, \bold{14}, 301--314.

Stock, J.H. and Watson, M.W. (2007). \emph{Introduction to Econometrics}, 2nd ed. Boston: Addison Wesley.
}

\seealso{\code{\link{StockWatson2007}}}

\examples{
data("CartelStability")
summary(CartelStability)
}

\keyword{datasets}
