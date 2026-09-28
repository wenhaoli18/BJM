# Build the prediction-to-infinity integration grid and marginal survival

Shared helper for `predictRisk` and `dynamicPredictionBio`: builds the
numerical-integration grid from `prediction_time` out to `upper_bound`,
and evaluates the marginal survival function `S(T)` over it.

## Usage

``` r
prepare_infinity_grid(
  data_predict_all,
  long_fit_all,
  survival_fit_all,
  prediction_time,
  upper_bound,
  bandcount2
)
```

## Value

A list with `predict.time.infinity`, `predict.time.infinity.1`, and
`S_T_all_infinity`.
