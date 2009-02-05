\name{ResumeNames}
\alias{ResumeNames}

\title{Are Emily and Greg More Employable Than Lakisha and Jamal?}

\description{
Cross-section data about resume, call-back and employer information for 4,870 fictitious
resumes.
}

\usage{data("ResumeNames")}

\format{
A data frame containing 4,870 observations on 27 variables.
\describe{
  \item{name}{factor indicating applicant's first name.}
  \item{gender}{factor indicating gender.}
  \item{ethnicity}{factor indicating ethnicity (i.e., Caucasian-sounding
    vs. African-American sounding first name).}
  \item{quality}{factor indicating quality of resume.}
  \item{call}{factor. Was the applicant called back?}
  \item{city}{factor indicating city: Boston or Chicago.}
  \item{jobs}{number of jobs listed on resume.}
  \item{experience}{number of years of work experience on the resume.}
  \item{honors}{factor. Did the resume mention some honors?}
  \item{volunteer}{factor. Did the resume mention some volunteering experience?}
  \item{military}{factor. Does the applicant have military experience?}
  \item{holes}{factor. Does the resume have some employment holes?}
  \item{school}{factor. Does the resume mention some work experience while at school?}
  \item{email}{factor. Was the e-mail address on the applicant's resume?}
  \item{computer}{factor. Does the resume mention some computer skills?}
  \item{special}{factor. Does the resume mention some special skills?}
  \item{college}{factor. Does the applicant have a college degree or more?}
  \item{minimum}{factor inidacting minimum experience requirement of the employer.}
  \item{equal}{factor. Is the employer EOE (equal opportunity employment)?}
  \item{wanted}{factor indicating type of position wanted by employer.}
  \item{requirements}{factor. Does the ad mention some requirement for the job?}
  \item{reqexp}{factor. Does the ad mention some experience requirement?}
  \item{reqcomm}{factor. Does the ad mention some communication skills requirement?}
  \item{reqeduc}{factor. Does the ad mention some educational requirement?}
  \item{reqcomp}{factor. Does the ad mention some computer skills requirement?}
  \item{reqorg}{factor. Does the ad mention some organizational skills requirement?}
  \item{industry}{factor indicating type of employer industry.}
}
}

\details{
Cross-section data about resume, call-back and employer information for 4,870 fictitious
resumes sent in response to employment advertisements in Chicago and Boston in 2001, 
in a randomized controlled experiment conducted by Bertrand and Mullainathan (2004).
The resumes contained information concerning the ethnicity of the applicant.
Because ethnicity is not typically included on a resume, resumes were differentiated on
the basis of so-called \dQuote{Caucasian sounding names} (such as Emily Walsh or Gregory Baker)
and \dQuote{African American sounding names} (such as Lakisha Washington or Jamal Jones).
A large collection of fictitious resumes were created and the pre-supposed
ethnicity (based on the sound of the name) was randomly assigned to each resume.
These resumes were sent to prospective employers to see which resumes generated a phone call
from the prospective employer.
}

\source{
Online complements to Stock and Watson (2007). 

\url{http://wps.aw.com/aw_stock_ie_2/}
}

\references{
Bertrand, M. and Mullainathan, S. (2004). Are Emily and Greg More Employable Than Lakisha and
Jamal? A Field Experiment on Labor Market Discrimination. 
\emph{American Economic Review}, \bold{94}, 991--1013.

Stock, J.H. and Watson, M.W. (2007). \emph{Introduction to Econometrics}, 2nd ed. Boston: Addison Wesley.
}

\seealso{\code{\link{StockWatson2007}}}

\examples{
data("ResumeNames")
summary(ResumeNames)
prop.table(xtabs(~ ethnicity + call, data = ResumeNames), 1)
}

\keyword{datasets}
