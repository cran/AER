\name{MurderRates}
\alias{MurderRates}

\title{Determinants of Murder Rates in the United States}

\description{
Cross-section data on states in 1950.
}

\usage{data("MurderRates")}

\format{
A data frame containing 44 observations on 8 variables.
\describe{
  \item{rate}{Murder rate per 100,000 (FBI estimate, 1950).}
  \item{convictions}{Number of convictions divided by number of murders in 1950.}
  \item{executions}{Average number of executions during 1946--1950 divided by convictions in 1950.}
  \item{time}{Median time served (in months) of convicted murderers released in 1951.}   
  \item{income}{Median family income in 1949 (in 1,000 USD).}  
  \item{lfp}{Labor force participation rate in 1950 (in percent).}
  \item{noncauc}{Proportion of population that is non-Caucasian in 1950.}
  \item{southern}{Factor indicating region.}
}
}

\source{
Maddala (2001), Table 8.4, p. 330
}

\references{
Maddala, G.S. (2001). \emph{Introduction to Econometrics}, 3rd ed. New York: John Wiley.

McManus, W.S. (1985). Estimates of the Deterrent Effect of Capital Punishment: 
The Importance of the Researcher's Prior Beliefs. \emph{Journal of Political Economy}, 
\bold{93}, 417--425.

Stokes, H. (2004). On the Advantage of Using Two or More Econometric Software Systems to Solve the Same Problem.
\emph{Journal of Economic and Social Measurement}, \bold{29}, 307--320.
}

\examples{
data("MurderRates")

## Maddala (2001, pp. 331)
fm_lm <- lm(rate ~ . + I(executions > 0), data = MurderRates)
summary(fm_lm)

model <- I(executions > 0) ~ time + income + noncauc + lfp + southern
fm_lpm <- lm(model, data = MurderRates)
summary(fm_lpm)

## Binomial models. Note: southern coefficient
fm_logit <- glm(model, data = MurderRates, family = binomial)
summary(fm_logit)

fm_logit2 <- glm(model, data = MurderRates, family = binomial,
  control = list(epsilon = 1e-15, maxit = 50, trace = FALSE))
summary(fm_logit2)

fm_probit <- glm(model, data = MurderRates, family = binomial(link = "probit"))
summary(fm_probit)

fm_probit2 <- glm(model, data = MurderRates , family = binomial(link = "probit"),
  control = list(epsilon = 1e-15, maxit = 50, trace = FALSE))
summary(fm_probit2)

## Explanation: quasi-complete separation
with(MurderRates, table(executions > 0, southern))
}

\keyword{datasets}
