# Point estimates, standard errors and complete-data df of one fit

Helper for
[`poolLongitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/poolLongitudinalSub.md).
For a continuous biomarker's
[`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) fit: the fixed
effects, their standard errors, and `lme()`'s own df per coefficient
(`fixDF$X`) – a between-subject covariate such as age has far fewer than
the within-subject time terms, so applying the intercept's df to every
coefficient would overstate the df of between-subject effects. For an
ordinal biomarker's
[`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html) fit: the
thresholds followed by the slopes, their Wald standard errors, and an
infinite complete-data df (`clmm()`'s inference is asymptotic).
Previously only `lme()` fits were handled, and pooling a fit with an
ordinal biomarker failed.

## Usage

``` r
marker_estimates(fit)
```

## Arguments

- fit:

  One biomarker's fitted marginal model.

## Value

A list with named vectors `estimate`, `std_error` and `dfcom`.
