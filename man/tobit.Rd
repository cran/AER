\name{tobit}
\alias{tobit}
\alias{print.tobit}
\alias{summary.tobit}
\alias{print.summary.tobit}
\alias{formula.tobit}
\alias{model.frame.tobit}
\alias{update.tobit}
\alias{waldtest.tobit}
\alias{lrtest.tobit}
\alias{linear.hypothesis.tobit}
\alias{logLik.survreg}
\alias{deviance.survreg}
\alias{fitted.survreg}

\title{Tobit Regression}

\description{
Fitting and testing tobit regression models for censored data.
}

\usage{
  tobit(formula, left = 0, right = Inf, dist = "gaussian", subset = NULL, data = list(), \dots)
}

\arguments{
  \item{formula}{a symbolic description of a regression model of type
    \code{y ~ x1 + x2 + \dots}.}
  \item{left}{left limit for the censored dependent variable \code{y}.
    If set to \code{-Inf}, \code{y} is assumed not be left-censored.}
  \item{right}{right limit for the censored dependent variable \code{y}.
    If set to \code{Inf}, the default, \code{y} is assumed not be right-censored.}
  \item{dist}{assumed distribution for the dependent variable \code{y}.
    This is passed to \code{\link[stats]{survreg}}, see the respective man page for
    more details.}
  \item{subset}{a specification of the rows to be used.}
  \item{data}{a data frame containing the variables in the model.}
  \item{\dots}{further arguments passed to \code{\link[stats]{survreg}}.}
}

\details{
  The function \code{tobit} is a convenience interface to \code{\link[stats]{survreg}}
  (for survival regression, including censored regression) setting different
  defaults and providing a more convenient interface for specification
  of the censoring information.
  
  The default is the classical tobit model (Tobin 1958, Greene 2003) assuming 
  a normal distribution for the dependent variable with left-censoring at 0.
  
  Technically, the formula of type \code{y ~ x1 + x2 + \dots} passed to \code{tobit}
  is simply transformed into a formula suitable for \code{\link[stats]{survreg}}: This means
  the dependent variable is first censored and then wrapped into a \code{\link[stats]{Surv}}
  object containing the censoring information which is subsequently passed to 
  \code{\link[stats]{survreg}}, e.g., \code{Surv(ifelse(y <= 0, 0, y), y > 0, type = "left") ~ x1 + x2 + \dots}
  for the default settings.
}

\references{
Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Upper Saddle River, NJ: Prentice Hall.

Tobin, J. (1958). Estimation of Relationships for Limited Dependent Variables.
  \emph{Econometrica}, \bold{26}, 24--36.
}

\value{
  An object of class \code{"tobit"} inheriting from class \code{"survreg"}.
}

\examples{
data("Affairs")

## from Table 22.4 in Greene (2003)
fm.tobit <- tobit(affairs ~ age + yearsmarried + religiousness + occupation + rating, data = Affairs)
fm.tobit2 <- tobit(affairs ~ age + yearsmarried + religiousness + occupation + rating,
  right = 4, data = Affairs)

summary(fm.tobit)
summary(fm.tobit2)
}

\keyword{regression}
