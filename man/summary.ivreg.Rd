\name{summary.ivreg}
\alias{summary.ivreg}
\alias{print.summary.ivreg}
\alias{vcov.ivreg}
\alias{bread.ivreg}
\alias{estfun.ivreg}
\alias{anova.ivreg}
\alias{hatvalues.ivreg}
\alias{predict.ivreg}
\alias{terms.ivreg}
\alias{model.matrix.ivreg}


\title{Methods for Instrumental-Variable Regression}
\description{
  Methods to standard generics for instrumental-variable regressions
  fitted by \code{\link{ivreg}}.
}

\usage{
\method{summary}{ivreg}(object, vcov. = NULL, df = NULL, diagnostics = FALSE, \dots)
\method{anova}{ivreg}(object, object2, test = "F", vcov = NULL, \dots)

\method{terms}{ivreg}(x, component = c("regressors", "instruments"), \dots)
\method{model.matrix}{ivreg}(object, component = c("projected", "regressors", "instruments"), \dots)
}

\arguments{
  \item{object, object2, x}{an object of class \code{"ivreg"} as fitted by
    \code{\link{ivreg}}.}
  \item{vcov., vcov}{a specification of the covariance matrix of the estimated
    coefficients. This can be specified as a matrix or as a function yielding a matrix
    when applied to the fitted model. If it is a function it is also employed in the two
    diagnostic F tests (if \code{diagnostics = TRUE} in the \code{summary()} method).}
  \item{df}{the degrees of freedom to be used. By default this is set to
    residual degrees of freedom for which a t or F test is computed. Alternatively,
    it can be set to \code{Inf} (or equivalently \code{0}) for which a z or Chi-squared
    test is computed.}
  \item{diagnostics}{logical. Should diagnostic tests for the instrumental-variable
    regression be carried out? These encompass an F test of the first stage regression
    for weak instruments, a Wu-Hausman test for endogeneity, and a Sargan
    test of overidentifying restrictions (only if there are more instruments than
    regressors).}
  \item{test}{character specifying whether to compute the large sample Chi-squared
    statistic (with asymptotic Chi-squared distribution) or the finite sample F statistic
    (with approximate F distribution).}
  \item{component}{character specifying for which component of the
    terms or model matrix should be extracted. \code{"projected"} gives the matrix of
    regressors projected on the image of the instruments.}
  \item{\dots}{currently not used.}
}

\details{
\code{\link{ivreg}} is the high-level interface to the work-horse function \code{\link{ivreg.fit}},
a set of standard methods (including \code{summary}, \code{vcov}, \code{anova},
\code{hatvalues}, \code{predict}, \code{terms}, \code{model.matrix}, \code{bread},
\code{estfun}) is available.
}

\seealso{\code{\link{ivreg}}, \code{\link[stats:lmfit]{lm.fit}}}

\examples{
## data
data("CigarettesSW")
CigarettesSW$rprice <- with(CigarettesSW, price/cpi)
CigarettesSW$rincome <- with(CigarettesSW, income/population/cpi)
CigarettesSW$tdiff <- with(CigarettesSW, (taxs - tax)/cpi)

## model 
fm <- ivreg(log(packs) ~ log(rprice) + log(rincome) | log(rincome) + tdiff + I(tax/cpi),
  data = CigarettesSW, subset = year == "1995")
summary(fm, vcov = sandwich, df = Inf, diagnostics = TRUE)

## ANOVA
fm2 <- ivreg(log(packs) ~ log(rprice) | tdiff, data = CigarettesSW, subset = year == "1995")
anova(fm, fm2, vcov = sandwich, test = "Chisq")
}

\keyword{regression}
