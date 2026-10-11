# Model-based AUC and Brier score for interval-censored outcomes

Every subject at risk counts, as a case with probability `p_case` and a
control with probability `p_control`
([`window_status_probs()`](https://wenhaoli18.github.io/BJM/reference/window_status_probs.md)):
the expected squared error for the Brier score, and case/control pairs
of different subjects weighted by \\p_i q_j\\ for the AUC (Yang et al.,
2026). These weights come from the model being evaluated, so the
measures lean optimistic when that model is misspecified.

## Usage

``` r
performance_metrics_interval_model(risk, p_case, p_control)
```

## Arguments

- risk:

  Predicted risks of the event type evaluated.

- p_case, p_control:

  Case and control probabilities per subject.

## Value

As
[`performance_metrics()`](https://wenhaoli18.github.io/BJM/reference/performance_metrics.md),
with `n_cases` the expected number of cases.
