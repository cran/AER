\name{ivreg.fit}
\alias{ivreg.fit}

\title{Fitting Instrumental-Variable Regressions}
\description{
  Fit instrumental-variable regression by two-stage least squares. This
  is equivalent to direct instrumental-variables estimation when the number of
  instruments is equal to the number of predictors.
}
\usage{
ivreg.fit(x, y, z, weights, offset, \dots)
}

\arguments{
  \item{x}{regressor matrix.}
  \item{y}{vector with dependent variable.}
  \item{z}{instruments matrix.}
  \item{weights}{an optional vector of weights to be used in the fitting process.}
  \item{offset}{an optional offset that can be used to specify an a priori known
    component to be included during fitting.}       
  \item{\dots}{further arguments passed to \code{\link[stats:lmfit]{lm.fit}} or
    \code{\link[stats:lmfit]{lm.wfit}}, respectively.}
}

\details{
\code{\link{ivreg}} is the high-level interface to the work-horse function \code{ivreg.fit},
a set of standard methods (including \code{summary}, \code{vcov}, \code{anova},
\code{hatvalues}, \code{predict}, \code{terms}, \code{model.matrix}, \code{bread},
\code{estfun}) is available and described on \code{\link{summary.ivreg}}.

\code{ivreg.fit} is a convenience interface to \code{\link[stats:lmfit]{lm.fit}} (or \code{\link[stats:lmfit]{lm.wfit}})
for first projecting \code{x} onto the image of \code{z} and the running 
a regression of \code{y} onto the projected \code{x}.
}

\value{
  \code{ivreg.fit} returns an unclassed list with the following components:
    \item{coefficients}{parameter estimates.}
    \item{residuals}{a vector of residuals.}
    \item{fitted.values}{a vector of predicted means.}
    \item{weights}{either the vector of weights used (if any) or \code{NULL} (if none).}
    \item{offset}{either the offset used (if any) or \code{NULL} (if none).}
    \item{estfun}{a matrix containing the empirical estimating functions.}
    \item{n}{number of observations.}
    \item{nobs}{number of observations with non-zero weights.}
    \item{rank}{the numeric rank of the fitted linear model.}
    \item{df.residual}{residual degrees of freedom for fitted model.}
    \item{cov.unscaled}{unscaled covariance matrix for the coefficients.}
    \item{sigma}{residual standard error.}
}

\seealso{\code{\link{ivreg}}, \code{\link[stats:lmfit]{lm.fit}}}

\examples{
## data
data("CigarettesSW")
CigarettesSW <- transform(CigarettesSW,
  rprice  = price/cpi,
  rincome = income/population/cpi,
  tdiff   = (taxs - tax)/cpi
)

## high-level interface
fm <- ivreg(log(packs) ~ log(rprice) + log(rincome) | log(rincome) + tdiff + I(tax/cpi),
  data = CigarettesSW, subset = year == "1995")

## low-level interface
y <- fm$y
x <- model.matrix(fm, component = "regressors")
z <- model.matrix(fm, component = "instruments")
ivreg.fit(x, y, z)$coefficients
}

\keyword{regression}
