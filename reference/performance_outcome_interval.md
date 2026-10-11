# Observed outcome per subject for an interval-censored fit

The interval-censored counterpart of
[`performance_outcome()`](https://wenhaoli18.github.io/BJM/reference/performance_outcome.md):
each subject's event is known only to lie in `(L, R]` (`R = Inf` if no
event was detected), with its type when detected.

## Usage

``` r
performance_outcome_interval(data, survival_fit_all, id_variable)
```

## Arguments

- data:

  The first `data_predict_all` data frame.

- survival_fit_all:

  An interval-censored `survivalSub.BJM` object.

- id_variable:

  The subject ID column.

## Value

A `data.frame` with columns `id`, `L`, `R`, `followup` (end of
follow-up: `R` after an event, `L` otherwise) and `cause` (1 or 2,
matching `risk_prob_1`/`risk_prob_2`; `NA` without competing risks or
without an event).
