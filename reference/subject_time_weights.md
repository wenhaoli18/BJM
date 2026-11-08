# Grid weights of a subject's event time (and event type)

\\\log \[f(Y \mid T, D) P(D \mid T) P(T \in \text{cell})\]\\ on a grid
of cells tiling \\(L, R\]\\, or \\(L, \infty)\\ for \\R = \infty\\ (the
equal-probability grid of
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)).
Without competing risks there is one column; with them one column per
event type in `types` – the known one, or both when `D` is `NA`. Shared
by
[`imputeEventTime()`](https://wenhaoli18.github.io/BJM/reference/imputeEventTime.md)
and
[`observedLongLogLik()`](https://wenhaoli18.github.io/BJM/reference/observedLongLogLik.md).

## Usage

``` r
subject_time_weights(
  data_j,
  use_y,
  long_fit_all,
  survival_fit_all,
  L,
  R,
  D,
  time_variable,
  survival_variable_all,
  survival_trans_function,
  n_grid,
  conditionalYT_fun = conditionalYT,
  conditionalYDT_fun = conditionalYDT,
  edges = NULL
)
```

## Arguments

- edges:

  Optional cell edges for a finite \\(L, R\]\\, in place of `n_grid`
  equal cells.

## Value

A list with `edges`, `mass` (the survival model's cell probabilities),
`log_w` (cells x types) and `types`.
