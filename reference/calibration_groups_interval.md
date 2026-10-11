# Observed risk per risk group for interval-censored outcomes

The interval-censored counterpart of
[`calibration_groups()`](https://wenhaoli18.github.io/BJM/reference/calibration_groups.md):
within each group of predicted risk, the observed risk of an event of
type `k` in \\(s, e\]\\ for a subject known to be event-free at their
last visit \\V\\ is \\(F_k(e) - F_k(s)) / (1 - F(V))\\ with \\F\\ from
[`npmle_interval()`](https://wenhaoli18.github.io/BJM/reference/npmle_interval.md),
averaged over the group's subjects – the same quantity
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
predicts. No confidence interval is given.

## Usage

``` r
calibration_groups_interval(risk, L, R, cause, V, k, s, end, n_groups)
```

## Arguments

- risk:

  Predicted risks of event type `k`.

- L, R, cause:

  Observed outcome (see
  [`performance_outcome_interval()`](https://wenhaoli18.github.io/BJM/reference/performance_outcome_interval.md)).

- V:

  Each subject's last visit up to the landmark.

- k, s, end:

  Event type, landmark and window end.

- n_groups:

  Number of groups.

## Value

As
[`calibration_groups()`](https://wenhaoli18.github.io/BJM/reference/calibration_groups.md).
