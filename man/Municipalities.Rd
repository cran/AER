\name{Municipalities}
\alias{Municipalities}

\title{Municipal Expenditure Data}

\description{
Panel data set for 265 Swedish municipalities covering 9 years (1979-1987).
}

\usage{data("Municipalities")}

\format{
A data frame containing 2,385 observations on 5 variables.
\describe{
  \item{municipality}{factor with ID number for municipality.}
  \item{year}{factor coding year.}
  \item{expenditures}{total expenditures.}
  \item{revenues}{total own-source revenues.}
  \item{grants}{intergovernmental grants received by the municipality.}
}
}

\details{
Total expenditures contains both capital and current expenditures. 

Expenditures, revenues, and grants are expressed in million SEK. The 
series are deflated and in per capita form. The implicit deflator is a 
municipality-specific price index obtained by dividing total local 
consumption expenditures at current prices by total local consumption 
expenditures at fixed (1985) prices. 

The data are gathered by Statistics Sweden and obtained from 
Financial Accounts for the Municipalities (Kommunernas Finanser). 
}

\source{
Journal of Applied Econometrics Data Archive. 

\url{http://www.econ.queensu.ca/jae/2000-v15.4/dahlberg-johansson/}
}

\references{
Dahlberg, M., and Johansson, E. (2000). An Examination of the Dynamic Behavior of Local Governments Using GMM Bootstrapping Methods. 
\emph{Journal of Applied Econometrics}, \bold{15}, 401--416.

Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Upper Saddle River, NJ: Prentice Hall.
}

\seealso{\code{\link{Greene2003}}}

\examples{
## Greene (2003), Table 18.2
data("Municipalities")
summary(Municipalities)
}

\keyword{datasets}
