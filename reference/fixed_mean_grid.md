# A biomarker's fixed-effects mean at every integration grid point

The longitudinal mean \\X(t)\beta\\ depends on the event time \\t\\ only
through the survival-time column and its `survival_variable_all`
transformations. This evaluates it for all rows of `data_i` (all
patients) at all grid points `l_i` with one
[`model.frame()`](https://rdrr.io/r/stats/model.frame.html)/[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
call on the rows stacked once per grid point, instead of one call per
patient and grid point, whose fixed overhead dominated the run time of
the prediction functions. The terms, factor levels and contrasts are
those of the fit, as in `predict.lme()`, so the design matrix is the one
the model was fit with, including data-dependent bases such as
[`poly()`](https://rdrr.io/r/stats/poly.html) and
[`splines::ns()`](https://rdrr.io/r/splines/ns.html). Transformed terms
such as `log(years)` are evaluated at each grid point. Rows are
processed in chunks of at most `max_rows` stacked rows.

## Usage

``` r
fixed_mean_grid(
  data_i,
  lfit_i,
  xlev_i,
  l_i,
  survival_variable,
  survival_variable_all,
  survival_trans_function,
  fixed_values = list(),
  beta = lfit_i$coefficients$fixed,
  max_rows = 2e+05
)
```

## Arguments

- data_i:

  Prediction rows for the biomarker.

- lfit_i:

  The biomarker's fitted model (`terms`, `contrasts`,
  `coefficients$fixed`).

- xlev_i:

  Its factor levels, or `NULL`.

- l_i:

  Grid of event times.

- fixed_values:

  Named list of further columns to set (e.g. the event type).

- beta:

  Fixed effects; by default those of `lfit_i`.

- max_rows:

  Largest number of stacked rows per chunk.

## Value

A `nrow(data_i)` x `length(l_i)` matrix; a row with a missing covariate
is `NA`.
