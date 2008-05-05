coeftest.multinom <- function(x, vcov. = NULL, df = NULL, ...)
{
  ## extract coefficients
  cc <- coef(x)
  est <- as.vector(t(cc))
  names(est) <- as.vector(t(outer(rownames(cc), colnames(cc), paste, sep = ":")))

  ## process vcov.
  if(is.null(vcov.)) vc <- vcov(x) else {
      if(is.function(vcov.)) vc <- vcov.(x)
        else vc <- vcov.
  }
  se <- sqrt(diag(vc))
  tval <- as.vector(est)/se

  ## process degrees of freedom  
  if(is.null(df)) df <- Inf

  if(is.finite(df) && df > 0) {
    pval <- 2 * pt(abs(tval), df = df, lower.tail = FALSE)
    cnames <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    mthd <- "t"
  } else {
    pval <- 2 * pnorm(abs(tval), lower.tail = FALSE)
    cnames <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    mthd <- "z"
  }
  rval <- cbind(est, se, tval, pval)
  colnames(rval) <- cnames
  class(rval) <- "coeftest"
  attr(rval, "method") <- paste(mthd, "test of coefficients")
  return(rval)
}

coeftest.polr <- function(x, vcov. = NULL, df = NULL, ...)
{
  ## extract coefficients
  est <- c(x$coefficients, x$zeta)

  ## process vcov.
  if(is.null(vcov.)) vc <- vcov(x) else {
      if(is.function(vcov.)) vc <- vcov.(x)
        else vc <- vcov.
  }
  se <- sqrt(diag(vc))
  tval <- as.vector(est)/se

  ## process degrees of freedom  
  if(is.null(df)) df <- Inf

  if(is.finite(df) && df > 0) {
    pval <- 2 * pt(abs(tval), df = df, lower.tail = FALSE)
    cnames <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    mthd <- "t"
  } else {
    pval <- 2 * pnorm(abs(tval), lower.tail = FALSE)
    cnames <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    mthd <- "z"
  }
  rval <- cbind(est, se, tval, pval)
  colnames(rval) <- cnames
  class(rval) <- "coeftest"
  attr(rval, "method") <- paste(mthd, "test of coefficients")
  return(rval)
}
