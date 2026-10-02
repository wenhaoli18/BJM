# IPCW time-dependent AUC and Brier score at one landmark

IPCW time-dependent AUC and Brier score at one landmark

## Usage

``` r
performance_metrics(risk, time, status, cause, k, s, horizon)
```

## Arguments

- risk:

  Predicted risks of event type `k` in the window.

- time, status, cause:

  Observed outcome of the same subjects (see
  [`performance_outcome()`](https://wenhaoli18.github.io/BJM/reference/performance_outcome.md));
  all have `time > s`.

- k:

  Event type evaluated (ignored when `cause` is all `NA`, i.e. without
  competing risks).

- s, horizon:

  Landmark and window length.

## Value

A list with `auc`, `brier`, `n_at_risk`, `n_cases`.
