# AER 1.2-15

* Updated reference output in some tests for CRAN checks.


# AER 1.2-14

* Updated reference output in some tests for CRAN checks.


# AER 1.2-13

* Fix `bread()` method for `tobit()` objects in case the latter already used a
  "robust" sandwich variance. In that case naive.var rather than var has
  to be used for the bread (reported by Daniel Klinenberg).


# AER 1.2-12

* Omitted some examples and updated corresponding output for CRAN checks.


# AER 1.2-11

* Rather than quitting from examples in the manual pages, when the required packages are
  not available, the examples now `stop()` in non-check settings.


# AER 1.2-10

* Improvements in `coeftest()` methods for `multinom` and `polr` objects to be compatible
  with recent extensions in `lmtest` >= 0.9-38.

* Fixed a bug in the tests for weak instruments in `ivdiag()`: Interactions like `x1:x2`
  and `x2:x1` were considered to be different variables although they are in fact equivalent.
  (Reported by Arne Henningsen and Ghislain Dossou.)

* Use fully-qualified calls to `survival::survreg()` and `survival::Surv()` within the
  `tobit()` function so that `AER::tobit()` can be called without attaching the package
  (reported by Bill Denney).

* Started replacing examples with repeated `with()` variable transformations with a
  single `transform()`.


# AER 1.2-9

* Some examples from the AER book lead to small numeric differences across different
  platforms. Hence, these have been excluded now from `tests/Ch-*.R` with an `IGNORE_RDIFF`
  comment for CRAN.

* The examples on the manual pages for the books `CameronTrivedi1998`, `Greene2003`,
  `WinkelmannBoes2009` are now excluded from testing to reduce the computational demands
  on CRAN. The corresponding analyses are mostly available on the manual pages for the
  underlying data sets, though.


# AER 1.2-8

* Starting from `survival` >= 3.1-6 the survival provides (or corrected) some standard
  S3 methods for `survreg` objects: `fitted()`, `nobs()`, `weights()`, `vcov()`. Previously,
  these were registered by AER" in order to work for "tobit" objects. For now,
  AER registers the methods for `tobit` objects. In the future, the methods will
  be dropped by AER altogether and inherited from survival (when AER depends on
  survival >= 3.1-6).


# AER 1.2-7

* Diagnostic tests `ivdiag()` for _weighted_ instrumental variables regression
  was incorrect as weights were erroneously ignored (detected and tested by
  Jonathan Siverskog).

* `vcov()` method for `survreg` objects in `survival` currently (version 2.44-1.1)
  provides no row/column names. Hence, a new `vcov()` method for `tobit` objects
  was added in `AER` mimicking the naming conventions from `survival:::summary.survreg`.
  Analogously, a `bread()` method was added for `tobit` objects that calls
  the `vcov()` method internally.

* `linearHypothesis()` method for `tobit()` objects no longer assumes that a scale
  parameter was estimated as part of the model. This enables the somewhat
  unusual but possible case of `summary(tobit(..., dist = "exponential"))`.

* In examples/demos/tests with calls to `set.seed(...)` the RNG version is now fixed
  to `suppressWarnings(RNGversion("3.5.0"))` to keep results exactly reproducible
  after fixing RNG problems in R 3.6.0.

* New errata item in `vignette("AER", package = "AER")`: The formula interface for
  `pgmm()` has changed. As of `plm` version 1.7-0, the function `dynformula()` is
  deprecated. Examples, tests, and demos have been adapted accordingly.

* In `tests/Ch-Intro.R` some `quantreg` computations are now wrapped into `try()`
  because the `summary()` method yields non-finite standard errors on some platforms.
  
* In `?CameronTrivedi1998` one version of the negbin-negbin hurdle model is now
  wrapped into `try()` because evaluating the log-likelihood at the given start
  values becomes too unstable on some platforms.


# AER 1.2-6

* Added `update()` method for `ivreg` objects that correctly handles the
  two-part right-hand side of the formula (suggested by Matthieu Stigler).

* Use `pdata.frame()` instead of `plm.data()` as preparation for modeling
  with `plm` and certain `systemfit` functions (requires `systemfit`
  1.1-20).

* The `model.matrix()` method for `ivreg` objects erroneously dropped
  the dimension of single-column projected regressor matrices. This
  propagated to the `estfun()` method and hence lead to problems with
  sandwich covariances (reported by Justus Winkelmann).

* In a bug fix the `survival` package changed the internal structure of
  `survreg()$y` starting from `survival` 2.42-7. This leads to incorrect
  `summary()` output for `tobit` objects which has been worked around
  now. Thanks to Terry Therneau for pointing out the problem and suggesting
  a fix.


# AER 1.2-5

* Support for aliased coefficients in `ivreg()` (suggested and tested by
  Liviu Andronic).

* New data sets `GoldSilver`, `MotorCycles2` and `MSCISwitzerland`, taken from
  Franses, van Dijk, Opschoor (2014): "Time Series Models for Business and
  Economic Forecasting", 2nd ed. For replication of the corresponding
  examples, several packages were added to the list of 'suggested' packages 
  (including `fGarch`, `forecast`, `longmemo`, `rugarch`, `vars`).

* Small improvements in `DESCRIPTION`/`Imports` and `NAMESPACE` for `R CMD check`.


# AER 1.2-4

* Reference output updated for recent versions of R.


# AER 1.2-3

* Package `splines` is loaded explicitly if needed (rather than assuming
  that it is loaded along with `survival`).

* Some URLs in the manual pages had been outdated and are updated (or omitted)
  now.


# AER 1.2-2

* Another bug fix in the new `summary(ivreg_object, diagnostics = TRUE)`. If
  sandwich standard errors (or other `vcov`) were used, the chi-squared form
  rather than the F form of the diagnostic Wald tests was computed and hence
  the p-values were incorrect.
  
* If there is more than one endogenous variable,
  `summary(ivreg_object, diagnostics = TRUE)`
  now reports separate tests of weak instruments for each endogenous variable.

  
# AER 1.2-1

* Bug fix in the new `summary(ivreg_object, diagnostics = TRUE)`. If there
  is more than one endogenous variable, the degrees of freedom (and hence
  the associated p-values) for the Sargan test were too large.
  
* The examples employing `rgl` for 3d visualizations (e.g., for the `SIC33`
  data) are not tested anymore in `R CMD check` (as `rgl` currently has some
  problems on CRAN's checks for OS X).


# AER 1.2-0

* The `summary()` method for `ivreg` now has a `diagnostics = FALSE` argument.
  If set to `TRUE`, three diagnostic tests are performed: an F test of
  the first stage regression for weak instruments, a Wu-Hausman test
  for endogeneity, and a Sargan test of overidentifying restrictions
  (only if there are more instruments than regressors).

* Added new data set `EquationCitations` provided by Fawcett & Higginson
  (2012, PNAS).

* Changes in Depends/Imports/Suggests due to new CRAN check requirements.
  In particular, the `Formula` package is now only imported (but not
  loaded into the search path).


# AER 1.1-9

* Recompressed data sets in package to reduce file storage requirements.

* `ivreg()` failed when used without instruments. Now fixed.

* The `summary()` for `ivreg` displayed the degrees of freedom of the
  overall Wald test incorrectly (although the p-value was computed
  correctly).

* Some technical changes for new R 2.14.0, e.g., adding `Authors@R` in
  `DESCRIPTION`, recompressing data, etc.


# AER 1.1-8

* The hat values for instrumental variables regressions are now
  computed in the `hatvalues()` method and not within `ivreg.fit()`
  to save computation time for large data sets.

* Added `nobs()` method for `survreg` objects (and thus `tobit` objects).
  Modified `ivreg` objects so that default `nobs()` methods works.

* Labeling in `coeftest()` method for `multinom` objects with
  binary responses has been fixed.

* Example 21.4 in `?Greene2003` now employs the scaled regressor
  `fincome/10000`.


# AER 1.1-7

* Adapted some example in `?Greene2003` in order to work both with
  the old and new `dynlm` package. `dynlm()` now provides convenient
  support for linear time trends via `dynlm(y ~ trend(y))` etc.    


# AER 1.1-6

* Adapted code/examples/tests to new `car` version which has deprecated
  `linear.hypothesis()` in favor of `linearHypothesis()`.


# AER 1.1-5

* `CPS1985` now has 534 observations (not 533 as in prior releases),
  the original observation 1 had been omitted inadvertently. See
  also the errata in `vignette("AER", package = "AER")`.

* Data and examples for Winkelmann and Boes (2009),
  "Analysis of Microdata" (2nd ed.) have been added. For details and
  extensive (but not quite complete) replication code see
  `help("WinkelmannBoes2009")` as well as `help("GSS7402")` and
  `help("GSOEP9402")`.

* As announced in the changes for version 1.1-0 of the `AER` package,
  the variable `earnings` has now been removed from the
  `PSID1976` (aka Mroz) data. In 1.1-0 it was renamed to `wage`
  to avoid confusion with other data sets.

* The `coeftest()` method for `polr` objects used to return wrong
  standard errors (and hence wrong z tests) for the last intercept.
  This was caused by an inconsistency between the `summary()` and
  `vcov()` methods for `polr` objects which has been improved in
  recent versions of the `MASS` package. The correct results are
  computed by `coeftest()` for `polr` objects computed with `MASS`
  version >= 7.3-6. See also the errata in
  `vignette("AER", package = "AER")`.

* The paper describing the various versions of the `Grunfeld` data
  has been accepted for publication in the German Economic Review.
  An updated version of the manuscript and associated replication
  files -- mostly based on `data("Grunfeld", package = "AER")` --
  is available from <https://www.zeileis.org/grunfeld/>.

* Added `lrtest()` method for `fitdistr` objects with intelligible
  model name (instead of the usual `formula` for formula-based models).


# AER 1.1-4

* `ivreg()` now uses the `Formula` package (>= 0.2-0) for processing
  of its model formulas. However, this only affects the internal
  computations, the user interface has remained unchanged.
  
* Numerous spelling improvements in the documentation (thanks to
  the new `aspell()` functionality in base R).


# AER 1.1-3

* Added `PSID7682` data set which contains the full Cornwell & Rupert
  (1988) panel data for the years 1976-1982. This should be used
  for estimation of the Hausman-Taylor model in Exercise 6 from
  Chapter 3 (instead of `PSID1982` which does not provide panel data
  but only the cross-section for 1982). See the errata and the
  manual page for more details.

* Fixed overall Wald test in `summary()` for `tobit` models
  with intercept only.
  

# AER 1.1-2

* New errata item in `vignette("AER", package = "AER")`: The comment
  regarding the output from the Johansen test (p. 169) is in error.
  The null hypothesis of no cointegration is not rejected at the
  10% level (only at 15% level).

* Enhancements of the `CigarettesSW` examples from Stock & Watson.

* Fixed overall Wald test in `summary()` for `tobit` models
  without intercept.
  
* Improved `rgl` code in the `SIC33` example.

* The variable `gender` in the `Parade2005` data set was
  wrong for observation 70. It is now `"male"` (not `"female"`). 
  

# AER 1.1-1

* A new improved version of the `plm` package is available
  from CRAN (version 1.1-1). This fixes a bug in the `summary()`
  of `plm` objects, see the vignette/errata for details.
  Furthermore, there is now a `vcovHC()` method for `panelmodel`
  objects: It gives equivalent results to `pvcovHC()` but
  is now the recommended user interface and hence used
  in the examples of some manual pages (see e.g. `?Fatalities`).


# AER 1.1-0

* Some variable names in the `PSID1976` (aka Mroz) data
  have been renamed: `earnings` is now called `wage`
  (to avoid confusion with other data sets), the previous
  variable `wage` is renamed as `repwage` (reported wage).
  Currently, `earnings` is kept; it will be removed in future
  releases.

* Documentation of the `Grunfeld` data has been enhanced and
  updated. Much more details are available in a recent
  technical report: Kleiber and Zeileis (2008), "The
  Grunfeld Data at 50", available from
  <https://doi.org/10.57938/c1925b32-0b27-490f-b593-42da650c196c>.
  
* Multinomial logit examples using Yves Croissant's `mlogit`
  package have been added for the `TravelMode` and `BankWages`
  data sets.

* Vignette/errata updated.
  

# AER 1.0-1

* Small changes for R 2.8.0.


# AER 1.0-0

* Official version accompanying the release of the
  book (contains all code from the book in demos
  and tests)

* See the new `vignette("AER", package = "AER")`
  for an overview of the package and a list of errata.


# AER 0.9-0

* Release of the version used for compiling the final
  version of the book for Springer.


# AER 0.2-0

* First CRAN release of the `AER` package.
