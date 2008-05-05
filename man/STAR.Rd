\name{STAR}
\alias{STAR}

\title{Project STAR: Student-Teacher Achievement Ratio}

\description{
The Project STAR public access data set, assessing the effect
of reducing class size on test scores in the early grades.
}

\usage{data("STAR")}

\format{
A data frame containing 11,598 observations on 47 variables.
\describe{
  \item{gender}{factor indicating student's gender.}
  \item{ethnicity}{factor indicating student's ethnicity with levels
    \code{"cauc"} (Caucasian), \code{"afam"} (African-American), \code{"asian"} (Asian), 
    \code{"hispanic"} (Hispanic), \code{"amindian"} (American-Indian) or \code{"other"}.}
  \item{birth}{student's birth quarter (of class \code{\link[zoo]{yearqtr}}).}
  \item{stark}{factor indicating the STAR class type in kindergarten:
    regular, small, or regular-with-aide. \code{NA} indicates that no STAR class was attended.}
  \item{star1}{factor indicating the STAR class type in 1st grade:
    regular, small, or regular-with-aide. \code{NA} indicates that no STAR class was attended.}
  \item{star2}{factor indicating the STAR class type in 2nd grade:
    regular, small, or regular-with-aide. \code{NA} indicates that no STAR class was attended.}
  \item{star3}{factor indicating the STAR class type in 3rd grade:
    regular, small, or regular-with-aide. \code{NA} indicates that no STAR class was attended.}
  \item{readk}{total reading scaled score in kindergarten.}
  \item{read1}{total reading scaled score in 1st grade.}
  \item{read2}{total reading scaled score in 2nd grade.}
  \item{read3}{total reading scaled score in 3rd grade.}
  \item{mathk}{total math scaled score in kindergarten.}
  \item{math1}{total math scaled score in 1st grade.}
  \item{math2}{total math scaled score in 2nd grade.}
  \item{math3}{total math scaled score in 3rd grade.}
  \item{lunchk}{factor indicating whether the student qualified for free lunch in kindergarten.}
  \item{lunch1}{factor indicating whether the student qualified for free lunch in 1st grade.}
  \item{lunch2}{factor indicating whether the student qualified for free lunch in 2nd grade.}
  \item{lunch3}{factor indicating whether the student qualified for free lunch in 3rd grade.}
  \item{schoolk}{factor indicating school type in kindergarten:
    \code{"inner-city"}, \code{"suburban"}, \code{"rural"} or \code{"urban"}.}
  \item{school1}{factor indicating school type in 1st grade:
    \code{"inner-city"}, \code{"suburban"}, \code{"rural"} or \code{"urban"}.}
  \item{school2}{factor indicating school type in 2nd grade:
    \code{"inner-city"}, \code{"suburban"}, \code{"rural"} or \code{"urban"}.}
  \item{school3}{factor indicating school type in 3rd grade:
    \code{"inner-city"}, \code{"suburban"}, \code{"rural"} or \code{"urban"}.}
  \item{degreek}{factor indicating gighest degree of kindergarten teacher:
    \code{"bachelor"}, \code{"master"}, \code{"specialist"}, or \code{"master+"}.}
  \item{degree1}{factor indicating gighest degree of 1st grade teacher:
    \code{"bachelor"}, \code{"master"}, \code{"specialist"}, or \code{"phd"}.}
  \item{degree2}{factor indicating gighest degree of 2nd grade teacher:
    \code{"bachelor"}, \code{"master"}, \code{"specialist"}, or \code{"phd"}.}
  \item{degree3}{factor indicating gighest degree of 3rd grade teacher:
    \code{"bachelor"}, \code{"master"}, \code{"specialist"}, or \code{"phd"}.}
  \item{ladderk}{factor indicating teacher's career ladder level in kindergarten: \code{"level1"},
    \code{"level2"}, \code{"level3"}, \code{"apprentice"}, \code{"probation"} or \code{"pending"}.}
  \item{ladder1}{factor indicating teacher's career ladder level in 1st grade: \code{"level1"},
    \code{"level2"}, \code{"level3"}, \code{"apprentice"}, \code{"probation"} or \code{"noladder"}.}
  \item{ladder2}{factor indicating teacher's career ladder level in 2nd grade: \code{"level1"},
    \code{"level2"}, \code{"level3"}, \code{"apprentice"}, \code{"probation"} or \code{"noladder"}.}
  \item{ladder3}{factor indicating teacher's career ladder level in 3rd grade: \code{"level1"},
    \code{"level2"}, \code{"level3"}, \code{"apprentice"}, \code{"probation"} or \code{"noladder"}.}
  \item{experiencek}{years of teacher's total teaching experience in kindergarten.}
  \item{experience1}{years of teacher's total teaching experience in 1st grade.}
  \item{experience2}{years of teacher's total teaching experience in 2nd grade.}
  \item{experience3}{years of teacher's total teaching experience in 3rd grade.}
  \item{tethnicityk}{factor indicating teacher's ethnicity in kindergarten with levels
    \code{"cauc"} (Caucasian) or \code{"afam"} (African-American).}
  \item{tethnicity1}{factor indicating teacher's ethnicity in 1st grade with levels
    \code{"cauc"} (Caucasian) or \code{"afam"} (African-American).}
  \item{tethnicity2}{factor indicating teacher's ethnicity in 2nd grade with levels
    \code{"cauc"} (Caucasian) or \code{"afam"} (African-American).}
  \item{tethnicity3}{factor indicating teacher's ethnicity in 3rd grade with levels
    \code{"cauc"} (Caucasian), \code{"afam"} (African-American), or \code{"asian"} (Asian).}
  \item{systemk}{factor indicating school system ID in kindergarten.}
  \item{system1}{factor indicating school system ID in 1st grade.}
  \item{system2}{factor indicating school system ID in 2nd grade.}
  \item{system3}{factor indicating school system ID in 3rd grade.}
  \item{schoolidk}{factor indicating school ID in kindergarten.}
  \item{schoolid1}{factor indicating school ID in 1st grade.}
  \item{schoolid2}{factor indicating school ID in 2nd grade.}
  \item{schoolid3}{factor indicating school ID in 3rd grade.}
}
}

\details{
Project STAR (Student/Teacher Achievement Ratio) was a four-year
longitudinal class-size study funded by the Tennessee General Assembly
and conducted in the late 1980s by the State Department of Education. Over 7,000 students
in 79 schools were randomly assigned into one of three interventions:
small class (13 to 17 students per teacher), regular class (22 to 25
students per teacher), and regular-with-aide class (22 to 25 students
with a full-time teacher's aide).  Classroom teachers were also randomly
assigned to the classes they would teach. The interventions were
initiated as the students entered school in kindergarten and continued
through third grade.

The Project STAR public access data set contains data on test scores,
treatment groups, and student and teacher characteristics for
the four years of the experiment, from academic year 1985--1986 to
academic year 1988--1989. The test score data analyzed in this
chapter are the sum of the scores on the math and reading portion of the
Stanford Achievement Test. 

Stock and Watson (2007) obtained the data set from the
Project STAR Web site at \url{http://www.heros-inc.org/star.htm}.

The data is provided in wide format. Reshaping it into long format
is illustrated below. Note that the levels of the \code{degree}, \code{ladder}
and \code{tethnicity} variables differ slightly between kindergarten
and higher grades.
}

\source{
Online complements to Stock and Watson (2007). 

\url{http://wps.aw.com/aw_stock_ie_2/}
}

\references{
Stock, J.H. and Watson, M.W. (2007). \emph{Introduction to Econometrics}, 2nd ed. Boston: Addison Wesley.
}

\seealso{\code{\link{StockWatson2007}}}

\examples{
data("STAR")

## Stock and Watson, p. 488
fmk <- lm(I(readk + mathk) ~ stark, data = STAR)
fm1 <- lm(I(read1 + math1) ~ star1, data = STAR)
fm2 <- lm(I(read2 + math2) ~ star2, data = STAR)
fm3 <- lm(I(read3 + math3) ~ star3, data = STAR)

coeftest(fm3, vcov = sandwich)
plot(I(read3 + math3) ~ star3, data = STAR)

## Stock and Watson, p. 489
fmke <- lm(I(readk + mathk) ~ stark + experiencek, data = STAR)
coeftest(fmke, vcov = sandwich)

## reshape data from wide into long format
## 1. variables and their levels
nam <- c("star", "read", "math", "lunch", "school", "degree", "ladder",
  "experience", "tethnicity", "system", "schoolid")
lev <- c("k", "1", "2", "3")
## 2. reshaping
star <- reshape(STAR, idvar = "id", ids = row.names(STAR),
  times = lev, timevar = "grade", direction = "long",
  varying = lapply(nam, function(x) paste(x, lev, sep = "")))
## 3. improve variable names and type
names(star)[5:15] <- nam
star$id <- factor(star$id)
star$grade <- factor(star$grade, levels = lev, labels = c("kindergarten", "1st", "2nd", "3rd"))
rm(nam, lev)

## fit a single model nested in grade (equivalent to fmk, fm1, fm2, fmk)
fm <- lm(I(read + math) ~ 0 + grade/star, data = star)
coeftest(fm, vcov = sandwich)

## visualization
library("lattice")
bwplot(I(read + math) ~ star | grade, data = star)
}

\keyword{datasets}
