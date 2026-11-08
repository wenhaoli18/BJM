# Fit the backward joint model with censored subjects in the longitudinal sub-model

**Experimental.** The usual right-censored fit estimates the
longitudinal sub-model \\f(Y \mid T)\\ on the subjects with an observed
event only. That is valid when censoring is independent of \\(T, Y)\\,
but \\f(Y \mid T)\\ is then learned only over the event times seen among
the events, which are the early ones when many subjects are censored.
Prediction needs \\f(Y \mid T)\\ for every \\T\\ after the prediction
time, so for a patient doing well it extrapolates the event-time effect
beyond the data. With a large censored fraction this shifts the
predictive distribution of \\T\\ towards early events: risks are
overestimated, increasingly so at later landmarks, and in short
follow-up the ranking of patients can be lost too.

`fitImputedBJM()` also uses the censored subjects, who carry the
information on long event times. It fills in their \\T \> C\\ by
stochastic EM, then by multiple imputation, as
[`fitIntervalBJM`](https://wenhaoli18.github.io/BJM/reference/fitIntervalBJM.md)
does for interval censoring:

1.  fit the survival sub-model
    ([`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md))
    once, and the longitudinal sub-model on the subjects with an event
    (the usual fit);

2.  draw each censored subject's \\T\\ from \\f(T \mid T \> C, Y, X)
    \propto f(Y \mid T, X) f(T \mid X)\\
    ([`imputeEventTime()`](https://wenhaoli18.github.io/BJM/reference/imputeEventTime.md));
    with competing risks, draw the event type \\D\\ jointly, from \\f(Y
    \mid T, D) P(D \mid T) f(T)\\;

3.  refit
    [`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
    on all subjects and repeat for `n_burnin` iterations;

4.  keep iterating, saving every `n_thin`-th fit until `n_imputations`
    are saved, and pool them for prediction (see
    [`poolLongitudinalSub`](https://wenhaoli18.github.io/BJM/reference/poolLongitudinalSub.md)).

The marginal survival model is the usual one (or the one passed in
`survival_fit_all`) and is fit once. With competing risks the event-type
model \\P(D \mid T)\\ is by default refit on every draw, on the subjects
with an event and the censored subjects with their drawn \\(T, D)\\, and
averaged over the imputations. The draws of \\T\\ beyond the last
follow-up rely on the survival model's tail: the extrapolated one by
default, or one borrowed from an external cohort with longer follow-up
([`survivalSubExternal`](https://wenhaoli18.github.io/BJM/reference/survivalSubExternal.md)).

Using the censored subjects makes the event-time effect matter where it
was never observed, so a term linear in \\T\\ can drift over the long
event times now in the fit. A bounded effect, e.g. `I(exp(-years / 5))`
in place of `years` in every term of `long_sub_fixed`, levels off for
long \\T\\; choose its scale by comparing
[`observedLongLogLik`](https://wenhaoli18.github.io/BJM/reference/observedLongLogLik.md)
across candidate fits. On `pbc3` (landmarks 1-3 years, 3-year window),
the usual fit overestimates the cumulative incidence of death by 8-32\\
cross-validation and this fit with `exp(-years / 5)` by -1 to 7\\ with
the same AUC.

## Usage

``` r
fitImputedBJM(
  data_survival,
  data_fit_all,
  form_marginal_surv,
  long_sub_fixed,
  long_sub_random,
  time_variable,
  form_conditional_cr = NULL,
  survival_variable_all = NULL,
  survival_trans_function = NULL,
  n_burnin = 5,
  n_imputations = 5,
  n_thin = 2,
  n_grid = 40,
  baseline = c("spline", "piecewise"),
  df = NULL,
  survival_fit_all = NULL,
  refit_event_type = TRUE,
  seed = NULL
)
```

## Arguments

- data_survival:

  One row per subject: event or censoring time, event indicator, the
  survival covariates and the subject id.

- data_fit_all:

  As in
  [`longitudinalSub`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md):
  one data frame per biomarker (or one shared data frame), long format,
  for **every** subject, with the event (or censoring) time in the
  column used in `long_sub_fixed`. The censored subjects' event-time
  column, event-type column and the columns named in
  `survival_variable_all` are (re)filled from the draws.

- form_marginal_surv, form_conditional_cr:

  As in
  [`survivalSub`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md):
  `Surv(time, status) ~ covariates` and the optional event-type model.
  The event-time variable of `form_marginal_surv` must be the one used
  in `long_sub_fixed`.

- long_sub_fixed, long_sub_random:

  As in
  [`longitudinalSub`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md).

- time_variable:

  Name of the visit-time column.

- survival_variable_all, survival_trans_function:

  Transformed event-time columns, as in
  [`predictRisk`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
  (`NULL` if none).

- n_burnin:

  Stochastic-EM iterations before fits are saved.

- n_imputations:

  Number of saved fits.

- n_thin:

  Iterations between saved fits.

- n_grid:

  Grid points per subject for the draws.

- baseline, df:

  Baseline hazard of the survival model (see
  [`survivalSub`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md)).

- survival_fit_all:

  Optional right-censored survival sub-model to use instead of fitting
  one from `form_marginal_surv` and `form_conditional_cr`, e.g. from
  [`survivalSubExternal`](https://wenhaoli18.github.io/BJM/reference/survivalSubExternal.md).
  `form_marginal_surv` is still used to read each subject's time and
  status in `data_survival`.

- refit_event_type:

  With competing risks, refit the event-type model on every draw (see
  Description); `FALSE` keeps the one fit on the subjects with an event.

- seed:

  Optional integer seed (the global RNG state is restored).

## Value

An object of class `"fitImputedBJM"`, a list with `survival_fit_all`
(with competing risks and `refit_event_type`, its event-type
coefficients averaged over the imputations; `survival_fit_all_list` has
one per imputation); `long_fit_all` (the imputations combined, for
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
and friends); `long_fit_all_list` (one fit per imputation); `pooled`
([`poolLongitudinalSub`](https://wenhaoli18.github.io/BJM/reference/poolLongitudinalSub.md)
of the fits, when `n_imputations > 1`); `imputed_T` (and, with competing
risks, `imputed_D`), censored subjects x imputations; `trace` (each
iteration's fixed effects, to check the chain has settled); and
`settings`. [`print()`](https://rdrr.io/r/base/print.html) shows a short
summary.

## See also

[`observedLongLogLik`](https://wenhaoli18.github.io/BJM/reference/observedLongLogLik.md)
to compare fits;
[`fitIntervalBJM`](https://wenhaoli18.github.io/BJM/reference/fitIntervalBJM.md)
for interval-censored event times.

## Examples

``` r
# \donttest{
data(pbc3)
pbc3_surv <- pbc3[!duplicated(pbc3$id), ]
fit <- fitImputedBJM(pbc3_surv, pbc3, Surv(years, status3) ~ age + sex,
                     long_sub_fixed = serBilir ~ year + age + sex +
                       I(exp(-years / 5)) + I(year * exp(-years / 5)),
                     long_sub_random = ~ year | id, time_variable = "year",
                     n_burnin = 2, n_imputations = 2, n_thin = 1, seed = 1)
fit
#> 
#> Backward joint model with imputed censored event times (fitImputedBJM)
#> 
#>   169 subjects with an event, 143 censored with T imputed; 2 imputation(s) after 2 burn-in iterations, keeping every iteration
#> 
#> Longitudinal sub-model, pooled over the imputations (Rubin's rules):
#> 
#>   [1] Outcome 1
#>                         estimate std_error     fmi
#> (Intercept)              0.41810  0.257200 0.10160
#> year                     0.03166  0.012570 0.02159
#> age                     -0.01279  0.004084 0.11400
#> sex                     -0.05711  0.129900 0.05040
#> I(exp(-years/5))         2.61600  0.155200 0.02854
#> I(year * exp(-years/5))  0.70130  0.052980 0.06498
#> 
#> Use x$long_fit_all and x$survival_fit_all for prediction; x$trace to check the chain has settled.

# compare with the usual fit on the subjects with an event
usual <- longitudinalSub(pbc3[pbc3$status3 == 1, ],
                         list(serBilir ~ year + age + sex + years + years * year),
                         list(~ year | id))
observedLongLogLik(fit, pbc3_surv, pbc3, time_variable = "year")
#> [1] -1356.182
observedLongLogLik(list(long_fit_all = usual, survival_fit_all = fit$survival_fit_all),
                   pbc3_surv, pbc3, time_variable = "year")
#> [1] -1573.894

history <- pbc3[pbc3$id == 2 & pbc3$year <= 3, ]
predictRisk(history, fit$long_fit_all, fit$survival_fit_all,
            prediction_time = 3, horizon = 3, time_variable = "year",
            survival_variable_all = NULL, survival_trans_function = NULL,
            bandcount1 = 10, bandcount2 = 10)
#> 
#> =================================================================
#>  Dynamic Prediction - Event Risk
#> -----------------------------------------------------------------
#>   Competing risks   : No
#>   Subjects          : 1
#> -----------------------------------------------------------------
#> 
#>  Subject Risk Prob
#>        2    0.1777
#> 
#> =================================================================
#> 
# }
```
