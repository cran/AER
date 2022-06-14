coeftest.multinom <- function(x, vcov. = NULL, df = NULL, ..., save = FALSE)
{
  ## extract coefficients
  est <- coef(x)
  if(!is.null(dim(est))) {
    est <- structure(as.vector(t(est)), 
      names = as.vector(t(outer(rownames(est), colnames(est), paste, sep = ":"))))
  }

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
  attr(rval, "df") <- df
  attr(rval, "logLik") <- logLik(x)
  if(save) attr(rval, "object") <- x
  return(rval)
}

coeftest.polr <- function(x, vcov. = NULL, df = NULL, ..., save = FALSE)
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
  attr(rval, "df") <- df
  attr(rval, "nobs") <- nobs(x)
  attr(rval, "logLik") <- logLik(x)
  if(save) attr(rval, "object") <- x
  return(rval)
}

lrtest.fitdistr <- function(object, ..., name = NULL)
{
  if(is.null(name)) name <- function(x) if(is.null(names(x$estimate))) {
    paste(round(x$estimate, digits = max(getOption("digits") - 3, 2)), collapse = ", ")
  } else {
    paste(names(x$estimate), "=", round(x$estimate, digits = max(getOption("digits") - 3, 2)), collapse = ", ")
  }
  lrtest.default(object, ..., name = name)
}

