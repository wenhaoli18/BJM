# Predicted and observed risk per risk group at one landmark

Predicted and observed risk per risk group at one landmark

## Usage

``` r
calibration_groups(risk, time, status, cause, k, end, n_groups)
```

## Arguments

- risk:

  Predicted risks of event type `k` in the window.

- time, status, cause:

  Observed outcome of the same subjects (see
  [`performance_outcome()`](https://wenhaoli18.github.io/BJM/reference/performance_outcome.md));
  all have `time` after the landmark.

- k:

  Event type evaluated (ignored without competing risks, i.e. when
  `cause` is all `NA`).

- end:

  End of the prediction window, `landmark + horizon`.

- n_groups:

  Number of groups (fewer if there are fewer subjects).

## Value

A `data.frame` with one row per group: `group`, `predicted` (mean
predicted risk), `observed`, `lower`, `upper`
(Kaplan–Meier/Aalen–Johansen estimate at `end` with 95\\
