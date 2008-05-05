\name{TechChange}
\alias{TechChange}
\title{Technological Change Data}

\description{
US time-series data, 1909--1949.  
}

\usage{data("TechChange")}

\format{
An annual multiple time series from 1909 to 1949 with 3 variables.
\describe{
  \item{output}{Output.}
  \item{clr}{Capital/labor ratio.}
  \item{technology}{Index of technology.}
}
}

\source{
Online complements to Greene (2003), Table F7.2.

\url{http://pages.stern.nyu.edu/~wgreene/Text/tables/tablelist5.htm}
}

\references{
Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Upper Saddle River, NJ: Prentice Hall.

Solow, R. (1957). Technical Change and the Aggregate Production Function. 
\emph{Review of Economics and Statistics}, \bold{39}, 312--320.
}

\seealso{\code{\link{Greene2003}}}

\examples{
data("TechChange")

## Greene (2003)
## Exercise 7.1
fm1 <- lm(I(output/technology) ~ log(clr), data = TechChange)
fm2 <- lm(I(output/technology) ~ I(1/clr), data = TechChange)
fm3 <- lm(log(output/technology) ~ log(clr), data = TechChange)
fm4 <- lm(log(output/technology) ~ I(1/clr), data = TechChange)

## Exercise 7.2 (a) and (c)
plot(I(output/technology) ~ clr, data = TechChange)
sctest(I(output/technology) ~ log(clr), data = TechChange, type = "Chow", point = c(1942, 1))
}

\keyword{datasets}
