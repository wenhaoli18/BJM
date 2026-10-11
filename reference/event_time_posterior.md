# Posterior distribution of the event time at a landmark

For the model-based performance measures: each patient's distribution of
\\T\\ (and event type) given their history up to the landmark \\s\\ and
\\T \> V\\, their last visit – the same quantity
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
integrates – tabulated on the cells of
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)'s
grids: \\(V, s\]\\
([`interval_gap_grid()`](https://wenhaoli18.github.io/BJM/reference/interval_gap_grid.md))
and \\(s, upper\]\\
([`prepare_infinity_grid()`](https://wenhaoli18.github.io/BJM/reference/prepare_infinity_grid.md)).

## Usage

``` r
event_time_posterior(
  data_predict_all,
  long_fit_all,
  survival_fit_all,
  s,
  time_variable,
  survival_variable_all,
  survival_trans_function,
  n_inf = 100,
  n_gap = 20
)
```

## Arguments

- data_predict_all:

  Landmark data (history up to `s`).

- long_fit_all, survival_fit_all:

  The fitted sub-models.

- s:

  The landmark.

- time_variable, survival_variable_all, survival_trans_function:

  As for
  [`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md).

- n_inf, n_gap:

  Grid cells after and before `s`.

## Value

A named list (patient ids) of lists with `edges` (cell boundaries, from
\\V\\) and `mass` (cells x event types, summing to 1).
