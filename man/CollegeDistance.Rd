\name{CollegeDistance}
\alias{CollegeDistance}

\title{College Distance Data}

\description{
Cross-section data from the High School and Beyond survey conducted by the
Department of Education in 1980, with a follow-up in 1986. The survey included
students from approximately 1,100 high schools.
}

\usage{data("CollegeDistance")}

\format{
A data frame containing 4,739 observations on 14 variables.
\describe{
  \item{gender}{factor indicating gender.}
  \item{ethnicity}{factor indicating ethincity (African-American, Hispanic or other).}
  \item{score}{base year composite test score. These are achievement tests given to high
    school seniors in the sample.}
  \item{fcollege}{factor. Is the father a college graduate?}
  \item{mcollege}{factor. Is the mother a college graduate?}
  \item{home}{factor. Does the family own their home?}
  \item{urban}{factor. Is the school in an urban area?}
  \item{unemp}{county unemployment rate in 1980.}
  \item{wage}{state hourly wage in manufacturing in 1980.}
  \item{distance}{distance from 4-year college (in 10 miles).}
  \item{tuition}{average state 4-year college tuition (in 1000 USD).}
  \item{education}{number of years of education.}
  \item{income}{factor. Is the family income above USD 25,000 per year?}
  \item{region}{factor indicating region (West or other).}
}
}

\details{
Rouse (1995) computed years of education by assigning 12 years to all
members of the senior class. Each additional year of secondary education counted as a
one year. Students with vocational degrees were assigned 13 years, AA degrees were
assigned 14 years, BA degrees were assigned 16 years, those with some graduate
education were assigned 17 years, and those with a graduate degree were assigned 18
years.

Stock and Watson (2007) provide separate data files for the students from
Western states and the remaining students. \code{CollegeDistance} includes
both data sets, subsets are easily obtained (see also examples).
}

\source{
Online complements to Stock and Watson (2007). 

\url{http://wps.aw.com/aw_stock_ie_2/}
}

\references{
Rouse, C.E. (1995). Democratization or Diversion? The Effect of Community
Colleges on Educational Attainment. \emph{Journal of Business \& Economic Statistics},
\bold{12}, 217--224.

Stock, J.H. and Watson, M.W. (2007). \emph{Introduction to Econometrics}, 2nd ed. Boston: Addison Wesley.
}

\seealso{\code{\link{StockWatson2007}}}

\examples{
## exclude students from Western states
data("CollegeDistance")
cd <- subset(CollegeDistance, region != "west")
summary(cd)
}

\keyword{datasets}
