\name{Electricity1970}
\alias{Electricity1970}

\title{Cost Function of Electricity Producers 1970}

\description{
Cross-section data, at the firm level, on electric power generation.  
}

\usage{data("Electricity1970")}

\format{
A data frame containing 158 cross-section observations on 9 variables.
\describe{
  \item{cost}{total cost.}
  \item{output}{total output.}
  \item{labor}{wage rate.}
  \item{laborshare}{cost share for labor.}
  \item{capital}{capital price index.}
  \item{capitalshare}{cost share for capital.}
  \item{fuel}{fuel price.}
  \item{fuelshare}{cost share for fuel.}
}
}

\details{
The data are from Christensen and Greene (1976) and pertain to the year 1970. 
However, the file contains some extra observations, the holding companies. 
Only the first 123 observations are needed to replicate Christensen and Greene (1976).
}

\source{
Online complements to Greene (2003), Table F5.2.

\url{http://pages.stern.nyu.edu/~wgreene/Text/tables/tablelist5.htm}
}

\references{
Christensen, L. and Greene, W.H. (1976). Economies of Scale in U.S. Electric Power Generation.
\emph{Journal of Political Economy}, \bold{84}, 655--676.

Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Upper Saddle River, NJ: Prentice Hall.
}

\seealso{\code{\link{Greene2003}}, \code{\link{Electricity1955}}}

\examples{
data("Electricity1970")

## Greene (2003), Ex. 5.6: a generalized Cobb-Douglas cost function
fm <- lm(log(cost/fuel) ~ log(output) + I(log(output)^2/2) + 
  log(capital/fuel) + log(labor/fuel), data=Electricity1970[1:123,])

}

\keyword{datasets}
