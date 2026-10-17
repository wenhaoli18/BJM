# Per-subject pieces of the copula E-step that do not depend on D

For every subject: the ordinal rows' fixed-effect linear predictor, the
thresholds bracketing each observed category and the random-effects
design of those rows; and the continuous rows' residuals from the fixed
effects, their random-effects design and residual precisions. Also the
positions `ro`/`rc` of the ordinal and continuous biomarkers' random
effects in `D`.

## Usage

``` r
copula_estep_setup(thetaLong, l, biomarker_type, thresholds)
```

## Arguments

- thetaLong:

  List with the stacked fixed effects `beta`, the starting `D` and the
  residual variances `sigma2` (1 for an ordinal biomarker).

- l:

  Per-subject design bookkeeping built by
  [`longitudinalSubCopula()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSubCopula.md);
  `l$yik[[m]]` holds an ordinal biomarker's category codes.

- biomarker_type, thresholds:

  As stored in the fit.

## Value

A list with `subjects` (one entry per subject), `ro` and `rc`.
