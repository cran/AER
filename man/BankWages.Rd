\name{BankWages} 
\alias{BankWages} 
\title{Bank Wages} 

\description{Wages of employees of a US bank.
} 

\usage{data("BankWages")}

\format{
A data frame containing 474 observations on 4 variables.
\describe{
  \item{job}{Ordered factor indicating job category, with levels \code{"custodial"},
    \code{"admin"}  and \code{"manage"}.}
  \item{education}{Education in years.}
  \item{gender}{Factor indicating gender.}
  \item{minority}{Factor. Is the employee member of a minority?}
}
}

\source{
Online complements to Heij, de Boer, Franses, Kloek, and van Dijk (2004). 

\url{http://www.oup.com/uk/booksites/content/0199268010/datasets/ch6/xr614bwa.asc}
}

\references{
Heij, C., de Boer, P.M.C., Franses, P.H., Kloek, T. and van Dijk, H.K. (2004). 
\emph{Econometric Methods with Applications in Business and Economics}. 
Oxford: Oxford University Press.
}

\examples{
data("BankWages")

## exploratory analysis of job ~ education
## (tables and spine plots, some education levels merged)
xtabs(~ education + job, data = BankWages)
edcat <- factor(BankWages$education)
levels(edcat)[3:10] <- rep(c("14-15", "16-18", "19-21"), c(2, 3, 3))
tab <- xtabs(~ edcat + job, data = BankWages)
prop.table(tab, 1)
spineplot(tab, off = 0)
plot(job ~ edcat, data = BankWages, off = 0)

## fit multinomial model for male employees
library("nnet")
fm_mnl <- multinom(job ~ education + minority, data = BankWages,
  subset = gender == "male", trace = FALSE)
summary(fm_mnl)
confint(fm_mnl)

## same with mlogit package
library("mlogit")
fm_mlogit <- mlogit(job ~ 1 | education + minority, data = BankWages, subset = gender == "male", shape = "wide", choice = "job", reflevel = "custodial")
summary(fm_mlogit)
}

\keyword{datasets}
