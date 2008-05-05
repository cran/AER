\name{MASchools} 
\alias{MASchools} 
\title{Massachusetts Test Score Data} 

\description{The dataset contains data on test performance, school 
characteristics and student demographic backgrounds for school districts 
in Massachusetts.
} 

\usage{data("MASchools")}

\format{
A data frame containing 220 observations on 16 variables.
\describe{
  \item{district}{character. District code.}
  \item{municipality}{character. Municipality name.}
  \item{expreg}{Expenditures per pupil, regular.}
  \item{expspecial}{Expenditures per pupil, special needs.}
  \item{expbil}{Expenditures per pupil, bilingual.}
  \item{expocc}{Expenditures per pupil, occupational.}
  \item{exptot}{Expenditures per pupil, total.}   
  \item{scratio}{Students per computer.}
  \item{special}{Special education students (per cent).}
  \item{lunch}{Percent qualifying for reduced-price lunch.}
  \item{stratio}{Student-teacher ratio.}
  \item{income}{Per capita income.}
  \item{score4}{4th grade score (math + English + science).}
  \item{score8}{8th grade score (math + English + science).}  
  \item{salary}{Average teacher salary.}
  \item{english}{Percent of English learners.}  
}
}

\details{The Massachusetts data are district-wide averages for 
public elementary school districts in 1998. The test score is taken 
from the Massachusetts Comprehensive Assessment System (MCAS) test, 
administered to all fourth graders in Massachusetts public schools 
in the spring of 1998. The test is sponsored by the Massachusetts 
Department of Education and is mandatory for all public schools.  
The data analyzed here are the overall total score, which is the sum 
of the scores on the English, Math, and Science portions of the 
test. Data on the student-teacher ratio, the percent of students 
receiving a subsidized lunch and on the percent of students still 
learning english are averages for each elementary school district 
for the 1997--1998 school year and were obtained from the 
Massachusetts department of education. Data on average district 
income are from the 1990 US Census.
}

\source{
Online complements to Stock and Watson (2007). 

\url{http://wps.aw.com/aw_stock_ie_2}
}

\references{
Stock, J. H. and Watson, M. W. (2007). \emph{Introduction to Econometrics}, 2nd ed. Boston: Addison Wesley.
}

\seealso{\code{\link{StockWatson2007}}, \code{\link{CASchools}}}

\examples{
## Massachusetts
data("MASchools")

## compare with California
data("CASchools")
CASchools$stratio <- with(CASchools, students/teachers)
CASchools$score4 <- with(CASchools, (math + read)/2)

## Stock and Watson, parts of Table 9.1, p. 330
vars <- c("score4", "stratio", "english", "lunch", "income")
cbind(
  CA_mean = sapply(CASchools[, vars], mean),
  CA_sd   = sapply(CASchools[, vars], sd),
  MA_mean = sapply(MASchools[, vars], mean),
  MA_sd   = sapply(MASchools[, vars], sd))

## Stock and Watson, Table 9.2, p. 332, col. (1)
fm1 <- lm(score4 ~ stratio, data = MASchools)
coeftest(fm1, vcov = vcovHC(fm1, type = "HC1"))

## More examples, notably the entire Table 9.2, can be found in:
## help("StockWatson2007")
}

\keyword{datasets}
