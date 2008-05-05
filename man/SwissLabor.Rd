\name{SwissLabor}
\alias{SwissLabor}
\title{Swiss Labor Market Participation Data}
\description{
Cross-section data originating from the health survey SOMIPOPS
for Switzerland in 1981.  
}
\usage{data("SwissLabor")}

\format{
A data frame containing 872 observations on 7 variables.
\describe{
  \item{participation}{Factor. Did the individual participate in the labor force?}
  \item{income}{Logarithm of nonlabor income.}
  \item{age}{Age in decades (years divided by 10).}
  \item{education}{Years of formal education.}
  \item{youngkids}{Number of young children (under 7 years of age).}
  \item{oldkids}{Number of older children (over 7 years of age).}
  \item{foreign}{Factor. Is the individual a foreigner (i.e., not Swiss)?}
}
}

\source{
Journal of Applied Econometrics Data Archive. 

\url{http://qed.econ.queensu.ca/jae/1996-v11.3/gerfin/}
}

\references{
Gerfin, M. (1996). Parametric and Semi-Parametric Estimation of the Binary Response
Model of Labour Market Participation. \emph{Journal of Applied Econometrics},
\bold{11}, 321--339.
}

\examples{
data("SwissLabor")

### Gerfin (1996), Table I.
fm_probit <- glm(participation ~ . + I(age^2), data = SwissLabor, family = binomial(link = "probit"))
summary(fm_probit)

### alternatively
fm_logit <- glm(participation ~ . + I(age^2), data = SwissLabor, family = binomial)
summary(fm_logit)
}

\keyword{datasets}
