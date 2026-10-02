# One row per subject with the observed outcome, for performancePlot()

One row per subject with the observed outcome, for performancePlot()

## Usage

``` r
performance_outcome(data, survival_fit_all, id_variable)
```

## Arguments

- data:

  The first `data_predict_all` data frame.

- survival_fit_all:

  A `survivalSub.BJM` object.

- id_variable:

  The subject ID column.

## Value

A `data.frame` with columns `id`, `time`, `status` (1 = event, 0 =
censored) and `cause` (1 or 2 for the event type matching
`risk_prob_1`/`risk_prob_2` of
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md);
`NA` without competing risks or when censored).
