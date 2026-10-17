# A biomarker's mean at every grid point, continuous or ordinal

For the copula density functions:
[`fixed_mean_grid()`](https://wenhaoli18.github.io/BJM/reference/fixed_mean_grid.md)
for a continuous biomarker; for an ordinal one, the linear predictor
\\X(t)\beta\\ of its `clmm()` fit (no intercept; the thresholds take its
place) from `long_sub_fixed[[i]]`, with all grid points stacked into one
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) call. An
ordinal formula with a term whose values depend on the other rows
([`poly()`](https://rdrr.io/r/stats/poly.html),
[`splines::ns()`](https://rdrr.io/r/splines/ns.html)/`bs()`,
[`factor()`](https://rdrr.io/r/base/factor.html),
[`scale()`](https://rdrr.io/r/base/scale.html),
[`cut()`](https://rdrr.io/r/base/cut.html), which
[`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
warns about) is evaluated one grid point at a time, as before, so that
its value does not change.

## Usage

``` r
copula_mean_grid(
  i,
  data_i,
  long_fit_all,
  l_i,
  survival_variable,
  survival_variable_all,
  survival_trans_function,
  fixed_values = list()
)
```

## Value

A `nrow(data_i)` x `length(l_i)` matrix.
