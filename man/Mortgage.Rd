\name{Mortgage}
\alias{Mortgage}
\title{Fixed versus Adjustable Mortgages}
\description{
Cross-section data about fixed versus adjustable mortgages for
78 households.  
}
\usage{data("Mortgage")}
\format{
A data frame containing 78 observations on 16 variables.
\describe{
\item{rate}{Factor with levels \code{"fixed"} and \code{"adjustable"}.}
\item{age}{Age of the borrower.}
\item{school}{Years of schooling for the borrower.}
\item{networth}{Net worth of the borrower.}
\item{interest}{Fixed interest rate.}
\item{points}{Ratio of points paid on adjustable to fixed rate mortgages.}
\item{maturities}{Ratio of maturities on adjustable to fixed rate mortgages.}
\item{years}{Years at the present address.}
\item{married}{Factor. Is the borrower married?}
\item{first}{Factor. Is the borrower a first-time home buyer?}
\item{selfemp}{Factor. Is the borrower self-employed?}
\item{tdiff}{The difference between the 10-year treasury rate less the 1-year treasury rate.}
\item{margin}{The margin on the adjustable rate mortgage.}
\item{coborrower}{Factor. Is there a co-borrower?}
\item{liability}{Short-term liabilities.}
\item{liquid}{Liquid assets.}
}
}

\source{
The data is from Baltagi (2002) and available at

\url{http://www.springeronline.com/sgw/cda/frontpage/0,10735,4-165-2-107420-0,00.html}
}

\references{
Baltagi, B.H. (2002). \emph{Econometrics}, 3rd ed. Berlin, Springer. 

Dhillon, U.S., Shilling, J.D. and Sirmans, C.F. (1987). Choosing Between Fixed and
Adjustable Rate Mortgages. \emph{Journal of Money, Credit and Banking}, \bold{19}, 260--267.

}

\seealso{\code{\link{Baltagi2002}}}

\examples{
data("Mortgage")
plot(rate ~ interest, data = Mortgage, breaks = fivenum(Mortgage$interest))
plot(rate ~ margin, data = Mortgage, breaks = fivenum(Mortgage$margin))
plot(rate ~ coborrower, data = Mortgage)
}

\keyword{datasets}
