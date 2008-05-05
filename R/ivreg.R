ivreg <- function(formula, instruments, data, subset, na.action, weights, offset,
  contrasts = NULL, model = TRUE, y = TRUE, x = FALSE, ...)
{
  ## set up model.frame() call  
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  
  ## extended formula processing, obtain:
  ##   ffx:     y ~ x1 + x2                (main regression)
  ##   ffz:       ~ z1 + z2 + z3           (instruments: can be NULL)
  ##   ff:      y ~ x1 + x2 | z1 + z2 + z3 (compact model formula)
  ##   formula: y ~ x1 + x2 + z1 + z2 + z3 (formula for setting up model frame)
  if(length(formula[[3]]) > 1 && identical(formula[[3]][[1]], as.name("|")))
  {
    ## regressors X and instruments Z specified by a single formula (with |)    
    if(!missing(instruments)) {
      warning("formula with '|' specified, 'instruments' ignored")
      cl$instruments <- NULL
    }
    ff <- formula
    formula[[3]][1] <- call("+")
    mf$formula <- formula
    ffx <- . ~ .
    ffz <- ~ .
    ffx[[2]] <- ff[[2]]
    ffx[[3]] <- ff[[3]][[2]]
    ffz[[3]] <- ff[[3]][[3]]
    ffz[[2]] <- NULL
    if(any(sapply(unlist(as.list(ffz[[2]])), function(x) identical(x, as.name("."))))) {
      ffz <- eval(parse(text = sprintf( paste("%s -", deparse(ffx[[2]])), deparse(ffz) )))
    }
  } else {
    if(missing(instruments)) {
      ## no instruments: plain OLS
      ffx <- ff <- formula
      ffz <- NULL
    } else {
      ## also support old "instruments" interface
      ffx <- formula
      if(length(instruments) > 2) instruments[[2]] <- NULL
      ffz <- instruments
      ff <- . ~ . | .
      ff[[2]] <- ffx[[2]]
      ff[[3]][[2]] <- ffx[[3]]
      ff[[3]][[3]] <- ffz[[2]]
      formula <- ff
      formula[[3]][1] <- call("+")
      mf$formula <- formula
      cl$instruments <- NULL
      cl$formula <- ff
    }
  }
  
  ## call model.frame()
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  ## extract response, terms, model matrices
  Y <- model.response(mf, "numeric")
  mt <- terms(formula, data = data)
  mtX <- terms(ffx, data = data)
  X <- model.matrix(mtX, mf, contrasts)
  if(is.null(ffz)) {
    mtZ <- NULL
    Z <- NULL
  } else {
    mtZ <- terms(ffz, data = data)
    mtZ <- terms(update(mtZ, ~ .), data = data)
    Z <- model.matrix(mtZ, mf, contrasts)
  }

  ## weights and offset
  weights <- model.weights(mf)
  offset <- model.offset(mf)
  if(is.null(offset)) offset <- 0
  if(length(offset) == 1) offset <- rep(offset, NROW(Y))
  offset <- as.vector(offset)

  ## call default interface
  rval <- ivreg.fit(X, Y, Z, weights, offset, ...)

  ## enhance information stored in fitted model object
  rval$call <- cl
  rval$formula <- ff
  rval$terms <- list(regressors = mtX, instruments = mtZ, full = mt)
  rval$na.action <- attr(mf, "na.action")
  rval$levels <- .getXlevels(mt, mf)
  rval$contrasts <- list(regressors = attr(X, "contrasts"), instruments = attr(Z, "contrasts"))
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- list(regressors = X, instruments = Z, projected = rval$x)
    else rval$x <- NULL
      
  class(rval) <- "ivreg"
  return(rval)
}

ivreg.fit <- function(x, y, z, weights, offset, ...)
{
  ## model dimensions
  n <- NROW(y)
  p <- ncol(x)
  
  ## defaults
  if(missing(z)) z <- NULL
  if(missing(weights)) weights <- NULL
  if(missing(offset)) offset <- rep(0, n)
  
  ## sanity checks
  stopifnot(n == nrow(x))
  if(!is.null(z)) stopifnot(n == nrow(z))
  if(!is.null(weights)) stopifnot(n == NROW(weights))
  stopifnot(n == NROW(offset))
  
  ## project regressors x on image of instruments z
  if(!is.null(z)) {
    auxreg <- if(is.null(weights)) lm.fit(z, x, ...) else lm.wfit(z, x, weights, ...)
    xz <- as.matrix(auxreg$fitted.values)
    pz <- z %*% chol2inv(auxreg$qr$qr) %*% t(z)
    colnames(xz) <- colnames(x)
  } else {
    xz <- x
    pz <- diag(NROW(x))
    colnames(pz) <- rownames(pz) <- rownames(x)
  }
  
  ## main regression
  fit <- if(is.null(weights)) lm.fit(xz, y, offset = offset, ...)
    else lm.wfit(xz, y, weights, offset = offset, ...)
 
  ## model fit information
  yhat <- drop(x %*% fit$coefficients)
  names(yhat) <- names(y)
  res <- y - yhat
  ucov <- chol2inv(fit$qr$qr[1:p, 1:p, drop = FALSE])
  colnames(ucov) <- rownames(ucov) <- names(fit$coefficients)
  rss <- if(is.null(weights)) sum(res^2) else sum(weights * res^2)
  hat <- diag(x %*% ucov %*% t(x) %*% pz)
  names(hat) <- rownames(x)

  rval <- list(
    coefficients = fit$coefficients,
    residuals = res,    
    fitted.values = yhat,
    weights = weights,
    offset = if(identical(offset, rep(0, n))) NULL else offset,
    n = n,
    rank = fit$rank,
    df.residual = fit$df.residual,
    cov.unscaled = ucov,
    sigma = sqrt(rss/fit$df.residual),
    hatvalues = hat,
    x = xz
  )
  
  return(rval)
}
   
vcov.ivreg <- function(object, ...)
  object$sigma^2 * object$cov.unscaled
    
bread.ivreg <- function (x, ...) 
    x$cov.unscaled * x$n

estfun.ivreg <- function (x, ...) 
{
    xmat <- model.matrix(x)
    wts <- weights(x)
    if(is.null(wts)) wts <- 1
    res <- residuals(x)
    rval <- as.vector(res) * wts * xmat
    attr(rval, "assign") <- NULL
    attr(rval, "contrasts") <- NULL
    return(rval)
}

hatvalues.ivreg <- function(model, ...)
  model$hatvalues

terms.ivreg <- function(x, component = c("regressors", "instruments"), ...)
  x$terms[[match.arg(component)]]

model.matrix.ivreg <- function(object, component = c("projected", "regressors", "instruments"), ...) {
  component <- match.arg(component)
  if(!is.null(object$x)) rval <- object$x[[component]]
    else if(!is.null(object$model)) {
      X <- model.matrix(object$terms$regressors, object$model, contrasts = object$contrasts$regressors)
      Z <- if(is.null(object$terms$instruments)) NULL
        else model.matrix(object$terms$instruments, object$model, contrasts = object$contrasts$instruments)
      w <- weights(object)
      XZ <- if(is.null(Z)) X
        else if(is.null(w)) lm.fit(Z, X)$fitted.values else lm.wfit(Z, X, w)$fitted.values
      rval <- switch(component,
        "regressors" = X,
	"instruments" = Z,
	"projected" = XZ)
    } else stop("not enough information in fitted model to return model.matrix")
  return(rval)
}

predict.ivreg <- function(object, newdata, na.action = na.pass, ...)
{
  if(missing(newdata)) fitted(object)
  else {
    mf <- model.frame(delete.response(object$terms$full), newdata,
      na.action = na.action, xlev = object$levels)
    X <- model.matrix(delete.response(object$terms$regressors), mf,
      contrasts = object$contrasts$regressors)
    drop(X %*% object$coefficients)
  }
}

print.ivreg <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  cat("Coefficients:\n")
  print.default(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
  cat("\n")
  invisible(x)
}

summary.ivreg <- function(object, vcov. = NULL, df = NULL, ...)
{
  ## weighted residuals
  res <- object$residuals
  y <- object$fitted.values + res
  n <- NROW(res)
  w <- object$weights
  if(is.null(w)) w <- rep(1, n)
  res <- res * sqrt(w)

  ## R-squared
  rss <- sum(res^2)
  if(attr(object$terms$regressors, "intercept")) {
    tss <- sum(w * (y - weighted.mean(y, w))^2)
    dfi <- 1    
  } else {
    tss <- sum(w * y^2)
    dfi <- 0
  }
  r.squared <- 1 - rss/tss
  adj.r.squared <- 1 - (1 - r.squared) * ((n - dfi)/object$df.residual)
  
  ## degrees of freedom (for z vs. t test)
  if(is.null(df)) df <- object$df.residual
  if(!is.finite(df)) df <- 0
  if(df > 0 & (df != object$df.residual)) {
    df <- object$df.residual
  }

  ## covariance matrix
  if(is.null(vcov.)) 
      vc <- vcov(object)
  else {
      if(is.function(vcov.)) vc <- vcov.(object)
        else vc <- vcov.
  }
  
  ## Wald test of each coefficient
  cf <- coeftest(object, vcov. = vc, df = df, ...)
  attr(cf, "method") <- NULL
  class(cf) <- "matrix"
  
  ## Wald test of all coefficients
  Rmat <- if(attr(object$terms$regressors, "intercept"))
    cbind(0, diag(length(coef(object))-1)) else diag(length(coef(object)))
  waldtest <- linear.hypothesis(object, Rmat, vcov. = vcov., test = ifelse(df > 0, "F", "Chisq"))
  waldtest <- c(waldtest[2,3], waldtest[2,4], -waldtest[2,2], if(df > 0) waldtest[1,1] else NULL)
  
  rval <- list(
    call = object$call,
    terms = object$terms,
    residuals = res,
    weights <- object$weights,
    coefficients = cf,
    sigma = object$sigma,    
    df = c(object$rank, if(df > 0) df else Inf, object$rank), ## aliasing not handled yet
    r.squared = r.squared,
    adj.r.squared = adj.r.squared,
    waldtest = waldtest,
    vcov = vc)
    
  class(rval) <- "summary.ivreg"
  return(rval)
}
 
print.summary.ivreg <- function(x, digits = max(3, getOption("digits") - 3), 
    signif.stars = getOption("show.signif.stars"), ...)
{
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat(if(!is.null(x$weights) && diff(range(x$weights))) "Weighted ", "Residuals:\n", sep = "")      
  if(NROW(x$residuals) > 5) {
      nam <- c("Min", "1Q", "Median", "3Q", "Max")
      rq <- if(length(dim(x$residuals)) == 2) 
	  structure(apply(t(x$residuals), 1, quantile), dimnames = list(nam, dimnames(x$residuals)[[2]]))
      else structure(quantile(x$residuals), names = nam)
      print(rq, digits = digits, ...)
  } else {
      print(x$residuals, digits = digits, ...)
  }

  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)

  cat("\nResidual standard error:", format(signif(x$sigma, digits)),
    "on", x$df[2], "degrees of freedom\n")

  cat("Multiple R-Squared:", formatC(x$r.squared, digits = digits))
  cat(",\tAdjusted R-squared:", formatC(x$adj.r.squared, digits = digits),
    "\nWald test:", formatC(x$waldtest[1], digits = digits),
    "on", x$waldtest[3], if(length(x$waldtest) > 3) c("and", x$waldtest[4]) else NULL,       
    "DF,  p-value:", format.pval(x$waldtest[2], digits = digits), "\n\n")
  invisible(x)
}
    
anova.ivreg <- function(object, object2, test = "F", vcov = NULL, ...)
{
  rval <- waldtest(object, object2, test = test, vcov = vcov)
  if(is.null(vcov)) {
    head <- attr(rval, "heading")
    head[1] <- "Analysis of Variance Table\n"
    rss <- sapply(list(object, object2), function(x) sum(residuals(x)^2))
    dss <- c(NA, -diff(rss))
    rval <- cbind(rval, cbind("RSS" = rss, "Sum of Sq" = dss))[,c(1, 5, 2, 6, 3:4)]
    attr(rval, "heading") <- head
    class(rval) <- c("anova", "data.frame")
  }
  return(rval)
}


## Falls #Instr. = #Regressoren, dann ist
##   b = (Z'X)^{-1} Z'y
## und loest die Schaetzgleichung
##   Z' (y - X beta) = 0
## Fuer
##   cov(y) = Omega
## haben wir dann
##   cov(b) = (Z'X)^{-1} Z' Omega Z (X'Z)^{-1}
##   
## Allgemein:  
##   b = (X' P_Z X)^{-1} X' P_Z y
## mit Schaetzgleichung
##   X' P_Z (y - X beta) = 0
## wobei P_Z der uebliche Projektor ist (Hat-Matrix bzgl Z) und
##   cov(b) = (X' P_Z X)^{-1} X' P_Z Omega P_Z X (X' P_Z X)^{-1}
## ist. M.a.W. meat ist X' P_Z Omega P_Z X, und bread ist (X' P_Z X)^{-1}
## 
## Siehe auch
##   http://www.stata.com/support/faqs/stat/2sls.html
