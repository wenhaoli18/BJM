# Build the prediction-to-infinity integration grid and marginal survival

Shared helper for `predictRisk` and `dynamicPredictionBio`: builds the
numerical-integration grid from `prediction_time` out to `upper_bound`,
and evaluates the marginal survival function `S(T)` over it.

The `bandcount2 + 1` intervals are spaced by probability, not by time:
each holds an equal share of the at-risk patients' (averaged)
model-based conditional survival mass beyond `prediction_time`. The
upper bound sits far out in the tail (see
[`integration_upper_bound()`](https://wenhaoli18.github.io/BJM/reference/integration_upper_bound.md)),
so equally spaced intervals would spend most of the grid where there is
almost no mass and converge more slowly as `bandcount2` grows. The first
interval starts exactly at `prediction_time`; the previous equally
spaced grid started half an interval earlier, counting mass from before
the prediction time. Each grid point is its interval's midpoint.

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

A list with `predict.time.infinity` (the grid points),
`predict.time.infinity.1` (the interval edges), and `S_T_all_infinity`.
