# Interval edges holding equal conditional survival mass

Helper for
[`prepare_infinity_grid()`](https://wenhaoli18.github.io/BJM/reference/prepare_infinity_grid.md):
`n_intervals + 1` edges from `prediction_time` to `upper_bound` such
that the average, over at-risk patients, of \\S(t \mid x) / S(s \mid
x)\\ drops by the same amount across every interval. Falls back to
equally spaced edges if no patient has a usable linear predictor.

## Usage

``` r
equal_mass_edges(
  data_predict_all,
  long_fit_all,
  survival_fit_all,
  prediction_time,
  upper_bound,
  n_intervals
)
```

## Value

A strictly increasing numeric vector of length `n_intervals + 1`.
