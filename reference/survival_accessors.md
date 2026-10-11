# Accessors for the marginal survival sub-model

Everything the prediction code needs from the marginal survival model of
[`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md)
– the name of the event-time variable, each patient's linear predictor,
the baseline cumulative hazard (per stratum), and the last time it is
tabulated at – goes through these helpers, so the same prediction code
serves both a right-censored fit (a
[`coxph`](https://rdrr.io/pkg/survival/man/coxph.html) model,
`survival_fit_all$coxph_fit`) and an interval-censored fit (a
piecewise-constant-hazard proportional hazards model,
`survival_fit_all$ic_fit`; see
[`icphFit()`](https://wenhaoli18.github.io/BJM/reference/icphFit.md)).
For a right-censored fit each accessor returns exactly what the code
used to compute from `coxph_fit` inline.

## Usage

``` r
is_interval_censored(survival_fit_all)

survival_time_variable(survival_fit_all)

survival_lp(survival_fit_all, newdata)

survival_cum_basehaz(survival_fit_all)

survival_patient_strata(survival_fit_all, newdata)

survival_max_followup(survival_fit_all)
```

## Arguments

- survival_fit_all:

  Output of
  [`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md).

- newdata:

  One row per patient.

## Functions

- `is_interval_censored()`: `TRUE` for an interval-censored fit.

- `survival_time_variable()`: Name of the event-time variable that the
  longitudinal sub-model conditions on (for an interval-censored fit,
  the `event_time` column given to
  [`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md)).

- `survival_lp()`: Uncentered linear predictor, `NA` for a patient with
  a missing covariate.

- `survival_cum_basehaz()`: Baseline cumulative hazard: a data frame
  with columns `hazard` and `time` (and `strata` for a stratified Cox
  model).

- `survival_patient_strata()`: Each patient's stratum (as character), or
  `NULL` when the model is not stratified.

- `survival_max_followup()`: Largest follow-up time in the data the
  survival model was fit on (for interval-censored data, the largest
  finite interval endpoint).
