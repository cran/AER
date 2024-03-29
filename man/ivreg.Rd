\name{ivreg}
\alias{ivreg}
\alias{print.ivreg}

\title{Instrumental-Variable Regression}
\description{
  Fit instrumental-variable regression by two-stage least squares. This
  is equivalent to direct instrumental-variables estimation when the number of
  instruments is equal to the number of predictors.
}
\usage{
ivreg(formula, instruments, data, subset, na.action, weights, offset,
  contrasts = NULL, model = TRUE, y = TRUE, x = FALSE, \dots)
}

\arguments{
  \item{formula, instruments}{formula specification(s) of the regression
    relationship and the instruments. Either \code{instruments} is missing and
    \code{formula} has three parts as in \code{y ~ x1 + x2 | z1 + z2 + z3}
    (recommended) or \code{formula} is \code{y ~ x1 + x2} and \code{instruments}
    is a one-sided formula \code{~ z1 + z2 + z3} (only for backward compatibility).}
  \item{data}{an optional data frame containing the variables in the model. 
    By default the variables are taken from the environment of the \code{formula}.}
  \item{subset}{an optional vector specifying a subset of observations to be used in 
    fitting the model.}
  \item{na.action}{a function that indicates what should happen when the 
    data contain \code{NA}s. The default is set by the \code{na.action} option.}
  \item{weights}{an optional vector of weights to be used in the fitting process.}
  \item{offset}{an optional offset that can be used to specify an a priori known
    component to be included during fitting.}       
  \item{contrasts}{an optional list. See the \code{contrasts.arg} of 
    \code{\link[stats:model.matrix]{model.matrix.default}}.}
  \item{model, x, y}{logicals.  If \code{TRUE} the corresponding components of
    the fit (the model frame, the model matrices , the response) are returned.}
  \item{\dots}{further arguments passed to \code{\link{ivreg.fit}}.}
}

\details{
\code{ivreg} is the high-level interface to the work-horse function \code{\link{ivreg.fit}},
a set of standard methods (including \code{print}, \code{summary}, \code{vcov}, \code{anova},
\code{hatvalues}, \code{predict}, \code{terms}, \code{model.matrix}, \code{bread},
\code{estfun}) is available and described on \code{\link{summary.ivreg}}.

Regressors and instruments for \code{ivreg} are most easily specified in a formula
with two parts on the right-hand side, e.g., \code{y ~ x1 + x2 | z1 + z2 + z3},
where \code{x1} and \code{x2} are the regressors and \code{z1},
\code{z2}, and \code{z3} are the instruments. Note that exogenous
regressors have to be included as instruments for themselves. For
example, if there is one exogenous regressor \code{ex} and one endogenous
regressor \code{en} with instrument \code{in}, the appropriate formula
would be \code{y ~ ex + en | ex + in}. Equivalently, this can be specified as
\code{y ~ ex + en | . - en + in}, i.e., by providing an update formula with a
\code{.} in the second part of the formula. The latter is typically more convenient,
if there is a large number of exogenous regressors.
}

\value{
  \code{ivreg} returns an object of class \code{"ivreg"}, with the following components:
    \item{coefficients}{parameter estimates.}
    \item{residuals}{a vector of residuals.}
    \item{fitted.values}{a vector of predicted means.}
    \item{weights}{either the vector of weights used (if any) or \code{NULL} (if none).}
    \item{offset}{either the offset used (if any) or \code{NULL} (if none).}
    \item{n}{number of observations.}
    \item{nobs}{number of observations with non-zero weights.}
    \item{rank}{the numeric rank of the fitted linear model.}
    \item{df.residual}{residual degrees of freedom for fitted model.}
    \item{cov.unscaled}{unscaled covariance matrix for the coefficients.}
    \item{sigma}{residual standard error.}
    \item{call}{the original function call.}
    \item{formula}{the model formula.}
    \item{terms}{a list with elements \code{"regressors"} and \code{"instruments"}
      containing the terms objects for the respective components.}
    \item{levels}{levels of the categorical regressors.}
    \item{contrasts}{the contrasts used for categorical regressors.}
    \item{model}{the full model frame (if \code{model = TRUE}).}
    \item{y}{the response vector (if \code{y = TRUE}).}
    \item{x}{a list with elements \code{"regressors"}, \code{"instruments"}, \code{"projected"},
      containing the model matrices from the respective components
      (if \code{x = TRUE}). \code{"projected"} is the matrix of regressors projected
      on the image of the instruments.}
}

\references{    
    Greene, W. H. (1993)
    \emph{Econometric Analysis}, 2nd ed., Macmillan.
}

\seealso{\code{\link{ivreg.fit}}, \code{\link[stats]{lm}}, \code{\link[stats:lmfit]{lm.fit}}}

\examples{
## data
data("CigarettesSW", package = "AER")
CigarettesSW <- transform(CigarettesSW,
  rprice  = price/cpi,
  rincome = income/population/cpi,
  tdiff   = (taxs - tax)/cpi
)

## model 
fm <- ivreg(log(packs) ~ log(rprice) + log(rincome) | log(rincome) + tdiff + I(tax/cpi),
  data = CigarettesSW, subset = year == "1995")
summary(fm)
summary(fm, vcov = sandwich, df = Inf, diagnostics = TRUE)

## ANOVA
fm2 <- ivreg(log(packs) ~ log(rprice) | tdiff, data = CigarettesSW, subset = year == "1995")
anova(fm, fm2)
}

\keyword{regression}
