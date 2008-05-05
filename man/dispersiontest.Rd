\name{dispersiontest}
\alias{dispersiontest}
\title{Dispersion Test}
\description{
 Tests the null hypothesis of equidispersion in Poisson GLMs against
 the alternative of overdispersion and/or underdispersion.
}
\usage{
dispersiontest(object, trafo = NULL, alternative = c("greater", "two.sided", "less"))
}
\arguments{
  \item{object}{a fitted Poisson GLM of class \code{"glm"} as fitted
    by \code{\link{glm}} with family \code{\link{poisson}}.}
 \item{trafo}{a specification of the alternative (see also details),
   can be numeric or a (positive) function or \code{NULL} (the default).}
 \item{alternative}{a character string specifying the alternative hypothesis:
   \code{"greater"} corresponds to overdispersion, \code{"less"} to
   underdispersion and \code{"two.sided"} to either one.}
}
\details{
The standard Poisson GLM models the (conditional) mean
\eqn{\mathsf{E}[y] = \mu}{E[y] = mu} which is assumed to be equal to the
variance \eqn{\mathsf{VAR}[y] = \mu}{VAR[y] = mu}. \code{dispersiontest}
assesses the hypothesis that this assumption holds (equidispersion) against
the alternative that the variance is of the form:
  \deqn{\mathsf{VAR}[y] \quad = \quad \mu \; + \; \alpha \cdot \mathrm{trafo}(\mu).}{VAR[y] = mu + alpha * trafo(mu).}
Overdispersion corresponds to \eqn{\alpha > 0}{alpha > 0} and underdispersion to
\eqn{\alpha < 0}{alpha < 0}. The coefficient \eqn{\alpha}{alpha} can be estimated
by an auxiliary OLS regression and tested with the corresponding t (or z) statistic
which is asymptotically standard normal under the null hypothesis.

Common specifications of the transformation function \eqn{\mathrm{trafo}}{trafo} are
\eqn{\mathrm{trafo}(\mu) = \mu^2}{trafo(mu) = mu^2} or \eqn{\mathrm{trafo}(\mu) = \mu}{trafo(mu) = mu}.
The former corresponds to a negative binomial (NB) model with quadratic variance function
(called NB2 by Cameron and Trivedi, 2005), the latter to a NB model with linear variance
function (called NB1 by Cameron and Trivedi, 2005) or quasi-Poisson model with dispersion 
parameter, i.e.,
  \deqn{\mathsf{VAR}[y] \quad = \quad (1 + \alpha) \cdot \mu = \mathrm{dispersion} \cdot \mu.}{VAR[y] = (1 + alpha) * mu = dispersion * mu.}
By default, for \code{trafo = NULL}, the latter dispersion formulation is used in
\code{dispersiontest}. Otherwise, if \code{trafo} is specified, the test is formulated
in terms of the parameter \eqn{\alpha}{alpha}. The transformation \code{trafo} can either
be specified as a function or an integer corresponding to the function \code{function(x) x^trafo},
such that \code{trafo = 1} and \code{trafo = 2} yield the linear and quadratic formulations
respectively.
}

\value{An object of class \code{"htest"}.}

\references{
Cameron, A.C. and Trivedi, P.K. (1990). Regression-based Tests for Overdispersion in the Poisson Model.  
\emph{Journal of Econometrics}, \bold{46}, 347--364.

Cameron, A.C. and Trivedi, P.K. (1998). \emph{Regression Analysis of Count Data}. 
Cambridge: Cambridge University Press.

Cameron, A.C. and Trivedi, P.K. (2005). \emph{Microeconometrics: Methods and Applications}.
Cambridge: Cambridge University Press.
}

\seealso{\code{\link{glm}}, \code{\link{poisson}}, \code{\link[MASS]{glm.nb}}}

\examples{
data("RecreationDemand")
rd <- glm(trips ~ ., data = RecreationDemand, family = poisson)

## linear specification (in terms of dispersion)
dispersiontest(rd)
## linear specification (in terms of alpha)
dispersiontest(rd, trafo = 1)
## quadratic specification (in terms of alpha)
dispersiontest(rd, trafo = 2)
dispersiontest(rd, trafo = function(x) x^2)

## further examples
data("DoctorVisits")
dv <- glm(visits ~ . + I(age^2), data = DoctorVisits, family = poisson)
dispersiontest(dv)

data("NMES1988")
nmes <- glm(visits ~ health + age + gender + married + income + insurance,
  data = NMES1988, family = poisson)
dispersiontest(nmes)
}
\keyword{htest}
