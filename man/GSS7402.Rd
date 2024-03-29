\name{GSS7402}
\alias{GSS7402}

\title{US General Social Survey 1974--2002}

\description{
Cross-section data for 9120 women taken from every fourth
year of the US General Social Survey between 1974 and 2002
to investigate the determinants of fertility.
}

\usage{data("GSS7402")}

\format{
A data frame containing 9120 observations on 10 variables. 
\describe{
  \item{kids}{Number of children. This is coded as a numerical variable
    but note that the value \code{8} actually encompasses 8 or more children.}
  \item{age}{Age of respondent.}
  \item{education}{Highest year of school completed.}
  \item{year}{GSS year for respondent.}
  \item{siblings}{Number of brothers and sisters.}
  \item{agefirstbirth}{Woman's age at birth of first child.}
  \item{ethnicity}{factor indicating ethnicity.
    Is the individual Caucasian (\code{"cauc"}) or not (\code{"other"})?}
  \item{city16}{factor. Did the respondent live in a city
    (with population > 50,000) at age 16?}
  \item{lowincome16}{factor. Was the income below average at age 16?}
  \item{immigrant}{factor. Was the respondent (or both parents)
    born abroad?}
}
}

\details{
This subset of the US General Social Survey (GSS) for every fourth year
between 1974 and 2002 has been selected by Winkelmann and Boes (2009)
to investigate the determinants of fertility. To do so they typically
restrict their empirical analysis to the women for which the completed
fertility is (assumed to be) known, employing the common cutoff of 40 years.
Both, the average number of children borne to a woman and the probability
of being childless, are of interest.
}

\source{
Online complements to Winkelmann and Boes (2009).
}

\references{
Winkelmann, R., and Boes, S. (2009). \emph{Analysis of Microdata}, 2nd ed. Berlin and Heidelberg: Springer-Verlag.
}

\seealso{\code{\link{WinkelmannBoes2009}}}

\examples{
\dontshow{ if(!requireNamespace("lattice") ||
              !requireNamespace("effects") ||
              !requireNamespace("MASS")) {
  if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
    stop("not all packages required for the example are installed")
  } else q() }}
## completed fertility subset
data("GSS7402", package = "AER")
gss40 <- subset(GSS7402, age >= 40)

## Chapter 1
## exploratory statistics
gss_kids <- prop.table(table(gss40$kids))
names(gss_kids)[9] <- "8+"

gss_zoo <- as.matrix(with(gss40, cbind(
  tapply(kids, year, mean),
  tapply(kids, year, function(x) mean(x <= 0)),
  tapply(education, year, mean))))
colnames(gss_zoo) <- c("Number of children",
  "Proportion childless", "Years of schooling")
gss_zoo <- zoo(gss_zoo, sort(unique(gss40$year)))

## visualizations instead of tables
barplot(gss_kids,
  xlab = "Number of children ever borne to women (age 40+)",
  ylab = "Relative frequencies")

library("lattice")
trellis.par.set(theme = canonical.theme(color = FALSE))
print(xyplot(gss_zoo[,3:1], type = "b", xlab = "Year"))


## Chapter 3, Example 3.14
## Table 3.1
gss40$nokids <- factor(gss40$kids <= 0, levels = c(FALSE, TRUE), labels = c("no", "yes"))
gss40$trend <- gss40$year - 1974
nokids_p1 <- glm(nokids ~ 1, data = gss40, family = binomial(link = "probit"))
nokids_p2 <- glm(nokids ~ trend, data = gss40, family = binomial(link = "probit"))
nokids_p3 <- glm(nokids ~ trend + education + ethnicity + siblings,
  data = gss40, family = binomial(link = "probit"))
lrtest(nokids_p1, nokids_p2, nokids_p3)


## Chapter 4, Figure 4.4
library("effects")
nokids_p3_ef <- effect("education", nokids_p3, xlevels = list(education = 0:20))
plot(nokids_p3_ef, rescale.axis = FALSE, ylim = c(0, 0.3))


## Chapter 8, Example 8.11
kids_pois <- glm(kids ~ education + trend + ethnicity + immigrant + lowincome16 + city16,
  data = gss40, family = poisson)
library("MASS")
kids_nb <- glm.nb(kids ~ education + trend + ethnicity + immigrant + lowincome16 + city16,
  data = gss40)
lrtest(kids_pois, kids_nb)


## More examples can be found in:
## help("WinkelmannBoes2009")
}

\keyword{datasets}
