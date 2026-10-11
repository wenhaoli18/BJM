# IPCW AUC and Brier score for interval-censored outcomes at one landmark

Uses only subjects whose status in \\(s, e\]\\ is certain: cases of type
`k` (last negative visit \\L \ge s\\, event detected at \\R \le e\\),
the other event type detected the same way (a non-case in the Brier
score), and controls (\\L \ge e\\). They are weighted by the inverse
probability of still being followed up – at \\R\\ for an event, at \\e\\
for a control – from the Kaplan–Meier estimate of the end of follow-up
among the subjects at risk at \\s\\, as in Yang et al. (2026),
Biometrical Journal 68:e70108. Subjects whose interval straddles \\s\\
or \\e\\ get weight 0.

## Usage

``` r
performance_metrics_interval_ipcw(risk, L, R, followup, cause, k, s, horizon)
```

## Arguments

- risk:

  Predicted risks of event type `k`.

- L, R, followup, cause:

  Observed outcome (see
  [`performance_outcome_interval()`](https://wenhaoli18.github.io/BJM/reference/performance_outcome_interval.md)).

- k, s, horizon:

  Event type, landmark and window length.

## Value

As
[`performance_metrics()`](https://wenhaoli18.github.io/BJM/reference/performance_metrics.md),
plus `n_uncertain`.
