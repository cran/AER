\name{ProgramEffectiveness}
\alias{ProgramEffectiveness}
\title{Program Effectiveness Data}

\description{
Data used to study the effectiveness of a program.
}

\usage{data("ProgramEffectiveness")}

\format{
A data frame containing 32 cross-section observations on 4 variables.
 \describe{
    \item{grade}{Factor with levels \code{"increase"} and \code{"decrease"}.}
    \item{average}{Grade-point average.}
    \item{testscore}{Test score on economics test.}
    \item{participation}{Factor. Did the individual participate in the program?}
    }
}

\details{
The data are taken form Spencer and Mazzeo (1980) who examined whether a new method 
of teaching economics significantly influenced performance in later economics courses.
}

\source{
Online complements to Greene (2003).

\url{http://pages.stern.nyu.edu/~wgreene/Text/tables/tablelist5.htm}
}

\references{
Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Upper Saddle River, NJ: Prentice Hall.

Spector, L. and Mazzeo, M. (1980). Probit Analysis and Economic Education. 
\emph{Journal of Economic Education}, \bold{11}, 37--44.
}

\seealso{\code{\link{Greene2003}}}

\examples{
data("ProgramEffectiveness")

## Greene (2003), Table 21.1, col. "Probit"
fm_probit <- glm(grade ~ average + testscore + participation,
  data = ProgramEffectiveness, family = binomial(link = "probit"))
summary(fm_probit)
}

\keyword{datasets}
