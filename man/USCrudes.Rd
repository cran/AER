\name{USCrudes}
\alias{USCrudes}
\title{US Crudes Data}

\description{
Cross-section data originating from 99 US oil field postings. 
}

\usage{data("USCrudes")}

\format{
A data frame containing 99 observations on 3 variables.
\describe{
  \item{price}{Crude prices (USD/barrel).}
  \item{gravity}{Gravity (degree API).}
  \item{sulphur}{Sulphur (in \%).}
}
}

\source{
The data is from Baltagi (2002).
}

\references{
Baltagi, B.H. (2002). \emph{Econometrics}, 3rd ed. Berlin, Springer. 

}

\seealso{\code{\link{Baltagi2002}}}

\examples{
\dontshow{ if(!requireNamespace("scatterplot3d")) {
  if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
    stop("not all packages required for the example are installed")
  } else q() }}
data("USCrudes")
plot(price ~ gravity, data = USCrudes)
plot(price ~ sulphur, data = USCrudes)
fm <- lm(price ~ sulphur + gravity, data = USCrudes)

## 3D Visualization
library("scatterplot3d")
s3d <- scatterplot3d(USCrudes[, 3:1], pch = 16)
s3d$plane3d(fm, lty.box = "solid", col = 4)
}

\keyword{datasets}
