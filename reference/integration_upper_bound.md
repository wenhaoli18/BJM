# Upper limit of the prediction-to-infinity integration grid

Shared helper for `predictRisk` and `dynamicPredictionBio`. The
denominator of a dynamic prediction integrates over every event time
after `prediction_time`, out to infinity; the grid has to stop
somewhere, and it should stop where the probability left beyond it is
negligible. So the upper limit is the earliest time by which every
at-risk patient's model-based conditional survival \\S(t \mid x) / S(s
\mid x)\\, \\s\\ = `prediction_time`, has dropped below `tail_prob`,
using the same (linearly extrapolated, per-stratum) baseline cumulative
hazard as
[`marginalT()`](https://wenhaoli18.github.io/BJM/reference/marginalT.md).
It never falls below `min_upper`.

This replaces twice the largest *observed survival time of the patients
being predicted*, which used each patient's own future outcome (not
available at `prediction_time`), stopped the integral early – and
inflated the risk – for a patient whose event came soon after
`prediction_time`, and failed when that time was missing.

Beyond the last time in the data
[`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md)
was fit on, the baseline hazard is an extrapolation; if it is so flat
that the tail criterion is not met by `max_multiple` times that last
time, the bound is capped there, with a warning if some patient's
conditional survival at the cap is still above `warn_prob`.

## Usage

``` r
integration_upper_bound(
  data_predict_all,
  long_fit_all,
  survival_fit_all,
  prediction_time,
  min_upper = prediction_time,
  tail_prob = 1e-04,
  max_multiple = 20,
  warn_prob = 0.01
)
```

## Arguments

- data_predict_all:

  At-risk prediction data (list of data frames).

- long_fit_all:

  Output of
  [`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
  (for the id variable).

- survival_fit_all:

  Output of
  [`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md).

- prediction_time:

  The prediction (landmark) time.

- min_upper:

  The bound is at least this (e.g. `prediction_time + horizon`, so the
  denominator grid covers the prediction window).

- tail_prob:

  Remaining conditional survival probability treated as negligible.

- max_multiple:

  Cap, as a multiple of the last training time.

- warn_prob:

  Warn when the cap leaves more than this conditional survival
  probability unintegrated for some patient.

## Value

A single number.
