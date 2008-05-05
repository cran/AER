\name{Journals}
\alias{Journals}
\title{Economics Journal Subscription Data}
\description{
Subscriptions to economics journals at US libraries, for the year 2000.
}
\usage{data("Journals")}

\format{
A data frame containing 180 observations on 10 variables.
\describe{
  \item{title}{Journal title.}
  \item{publisher}{factor with publisher name.}
  \item{society}{factor. Is the journal published by a scholarly society?}
  \item{price}{Library subscription price.}
  \item{pages}{Number of pages.}
  \item{charpp}{Characters per page.}
  \item{citations}{Total number of citations.}
  \item{foundingyear}{Year journal was founded.}
  \item{subs}{Number of library subscriptions.}
  \item{field}{factor with field description.}
}
}

\details{
Data on 180 economic journals, collected in particular for analyzing journal
pricing. See also \url{http://www.econ.ucsb.edu/~tedb/Journals/jpricing.html}
for general information on this topic as well as a more up-to-date version of
the data set. This version is taken from Stock and Watson (2007).

The data as obtained from \url{http://wps.aw.com/aw_stock_ie_2}
contained two journals with title \dQuote{World Development}. One of these (observation 80)
seemed to be an error and was changed to \dQuote{The World Economy}.
}

\source{
Online complements to Stock and Watson (2007). 

\url{http://wps.aw.com/aw_stock_ie_2/0,12040,3332253-,00.html}
}

\references{
Bergstrom, T. (2001). Free Labor for Costly Journals? \emph{Journal of Economic Perspectives}, 15, 183--198.

Stock, J.H. and Watson, M.W. (2007). \emph{Introduction to Econometrics}, 2nd ed. Boston: Addison Wesley.
}

\seealso{\code{\link{StockWatson2007}}}

\examples{
## data and transformed variables
data("Journals")
journals <- Journals[, c("subs", "price")]
journals$citeprice <- Journals$price/Journals$citations
journals$age <- 2000 - Journals$foundingyear
journals$chars <- Journals$charpp*Journals$pages/10^6

## Stock and Watson (2007)
## Figure 8.9 (a) and (b)
plot(subs ~ citeprice, data = journals, pch = 19)
plot(log(subs) ~ log(citeprice), data = journals, pch = 19)
fm1 <- lm(log(subs) ~ log(citeprice), data = journals)
abline(fm1)

## Table 8.2, use HC1 for comparability with Stata 
fm2 <- lm(subs ~ citeprice + age + chars, data = log(journals))
fm3 <- lm(subs ~ citeprice + I(citeprice^2) + I(citeprice^3) +
  age + I(age * citeprice) + chars, data = log(journals))
fm4 <- lm(subs ~ citeprice + age + I(age * citeprice) + chars, data = log(journals))
coeftest(fm1, vcov = vcovHC(fm1, type = "HC1"))
coeftest(fm2, vcov = vcovHC(fm2, type = "HC1"))
coeftest(fm3, vcov = vcovHC(fm3, type = "HC1"))
coeftest(fm4, vcov = vcovHC(fm4, type = "HC1"))
waldtest(fm3, fm4, vcov = vcovHC(fm3, type = "HC1"))

## changes with respect to age
## Nyblom-Hansen test
scus <- gefp(subs ~ citeprice, data = log(journals), fit = lm, order.by = ~ age)
plot(scus, functional = meanL2BB)
## estimate breakpoint(s)
journals <- journals[order(journals$age),]
bp <- breakpoints(subs ~ citeprice, data = log(journals), h = 20)
plot(bp)
bp.age <- journals$age[bp$breakpoints]
## visualization
plot(subs ~ citeprice, data = log(journals), pch = 19, col = (age > log(bp.age)) + 1)
abline(coef(bp)[1,], col = 1)
abline(coef(bp)[2,], col = 2)
legend("bottomleft", legend = c("age > 18", "age < 18"), lty = 1, col = 2:1, bty = "n")
}

\keyword{datasets}
