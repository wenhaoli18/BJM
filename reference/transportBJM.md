# Transport a fitted backward joint model to a new population

**Experimental.** Adapts a survival sub-model fitted in one population
(the *source*, e.g. the hospital the model was developed in) to a
*target* population (e.g. a new hospital) whose event rates differ,
while keeping the longitudinal sub-model unchanged.

The backward joint model factors the joint distribution as \\f(T, D \mid
X)\\ f(Y \mid T, D, X)\\: the event-time model and the biomarker
trajectory *given* the event time. Between populations it is often the
first factor that changes most – a different case mix, standard of care,
or referral pattern moves the event rates – while how a biomarker
evolves in the years before an event is closer to biology and more
stable. This is the *label shift* assumption of machine learning, and
the backward factorization is exactly the one under which it is easy to
act on: only `survival_fit_all` needs to change, and the result plugs
into
[`predictRisk`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md),
[`predictLongitudinal`](https://wenhaoli18.github.io/BJM/reference/predictLongitudinal.md),
[`simulateTrajectory`](https://wenhaoli18.github.io/BJM/reference/simulateTrajectory.md),
[`performancePlot`](https://wenhaoli18.github.io/BJM/reference/performancePlot.md)
and
[`calibrationPlot`](https://wenhaoli18.github.io/BJM/reference/calibrationPlot.md)
with the source `long_fit_all`.

Three methods are available, by how much is known about the target:

- `"baseline"`:

  Target outcomes are available (`data_target`, one row per subject).
  The hazard ratios are kept from the source and only the baseline
  cumulative hazard is re-estimated on the target (Breslow estimator
  with the source linear predictor as an offset). Needs few target
  events.

- `"recalibrate"`:

  As `"baseline"`, but the source linear predictor is also rescaled by a
  calibration slope \\\gamma\\, estimated by a Cox model of the target
  outcome on the source linear predictor (\\\gamma \< 1\\: the source
  hazard ratios are too extreme for the target).

- `"unlabeled"`:

  No target outcomes at all, only biomarker histories
  (`data_predict_all`). The target baseline cumulative hazard is taken
  to be the source one times \\e^{\alpha}\\, and \\\alpha\\ is estimated
  by maximizing the marginal likelihood of the target histories,
  \$\$\prod_i \int\_{s_i}^{\infty} f(Y_i \mid T)\\ f\_\alpha(T \mid
  X_i)\\ dT,\$\$ where \\s_i\\ is the time up to which patient \\i\\ is
  known to be event-free (by default their last measurement). Because
  \\f(Y \mid T)\\ is held fixed, the shape of the target trajectories –
  how close they look to the trajectories of patients about to have an
  event – tells how far away the events are. This is the same integral
  [`predictRisk`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
  uses as its denominator, and it is only possible because the biomarker
  model conditions on the event time.

## Usage

``` r
transportBJM(
  survival_fit_all,
  method = c("baseline", "recalibrate", "unlabeled"),
  data_target = NULL,
  update_event_type = FALSE,
  data_predict_all = NULL,
  long_fit_all = NULL,
  time_variable = NULL,
  survival_variable_all = NULL,
  survival_trans_function = NULL,
  landmark = NULL,
  data_source = NULL,
  cap_event_time = TRUE,
  n_grid = 200,
  upper = NULL,
  alpha_range = c(-4, 4)
)
```

## Arguments

- survival_fit_all:

  Output of
  [`survivalSub`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md)
  fit on the source data (right-censored, not stratified).

- method:

  `"baseline"`, `"recalibrate"` or `"unlabeled"`; see Description.

- data_target:

  For `"baseline"` and `"recalibrate"`: one row per target subject, with
  the outcome and covariates used in
  `survival_fit_all$form_marginal_surv` (and, with
  `update_event_type = TRUE`, `form_conditional_cr`).

- update_event_type:

  For `"baseline"`/`"recalibrate"` with competing risks: also refit the
  event-type model on the target (default `FALSE`, keep the source one).

- data_predict_all, long_fit_all, time_variable, survival_variable_all,
  survival_trans_function:

  For `"unlabeled"`: the target biomarker histories and the source
  longitudinal sub-model, as for
  [`predictRisk`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md).
  The event-time column may be absent or `NA`.

- landmark:

  For `"unlabeled"`: the time up to which each target patient is known
  to be event-free, a single number or a vector named by patient id.
  Defaults to each patient's last measurement time.

- data_source:

  For `"unlabeled"`: optionally, biomarker histories from the source
  population (in the format of `data_predict_all`). The shift is then
  estimated on them too, where it should be 0, and subtracted: whatever
  misfit of the model biases the estimate in the same way in both
  populations cancels. Recommended.

- cap_event_time:

  For `"unlabeled"`: if `TRUE` (default), an event time after the source
  follow-up is given the biomarker trajectory of an event at the end of
  follow-up, as in
  [`simulateTrajectory`](https://wenhaoli18.github.io/BJM/reference/simulateTrajectory.md),
  instead of extrapolating the longitudinal sub-model.

- n_grid, upper:

  For `"unlabeled"`: the number of grid intervals on `(0, upper]` used
  for the integral, and the upper end of the grid (default three times
  the source follow-up); the probability of an event beyond `upper` is
  placed at `upper`.

- alpha_range:

  For `"unlabeled"`: the search interval for \\\alpha\\, the log ratio
  of target to source baseline hazard.

## Value

An object of class `c("transportBJM", "survivalSub.BJM")`:
`survival_fit_all` with an extra element `transport` (the method, the
target baseline cumulative hazard, the calibration slope, and `alpha`
with its standard error for `"unlabeled"`). Use it in place of
`survival_fit_all` anywhere in the package.

## Details

Only the marginal event-time model is transported; the longitudinal
sub-model, and by default the event-type model, are the source ones. If
the biomarker trajectories before an event themselves differ between the
populations (different assays, measurement schedules that depend on
health), label shift does not hold and `"unlabeled"` is biased.
`"unlabeled"` identifies a single proportional shift; it carries no
information from patients without any measurement, and little from
patients with only one.

## Examples

``` r
data(pbc3)
first <- pbc3[!duplicated(pbc3$id), ]
odd <- seq_len(nrow(first)) %% 2 == 1
source <- first[!odd, ]
target <- first[odd, ]
fit_src <- survivalSub(source, Surv(years, status3) ~ age + sex, NULL)
fit_tgt <- transportBJM(fit_src, "baseline", data_target = target)
fit_tgt
#> Backward joint model survival sub-model, transported to a target population
#>   Method                    : baseline
#>   Target subjects used      : 156
#>   Target events             : 84
#> 
#> Source fit:
#> 
#> Call:
#> survivalSub(form_marginal_surv = Surv(years, status3) ~ age + sex)
#> 
#> Data Descriptives:
#>   Number of subjects        : 156
#>   Number of events          : 85
#> 
#> =================================================================
#>  Marginal Survival Sub-model  [Cox PH]
#> -----------------------------------------------------------------
#>  Formula: Surv(years, status3) ~ age + sex
#> <environment: 0x61a5cb3c2548>
#> 
#>         Coef exp(Coef)       SE      z p-value  
#> age  0.01835   1.01852  0.01050  1.748  0.0804 .
#> sex -0.42317   0.65497  0.28575 -1.481  0.1386  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>   n = 156,  events = 85
#>   Concordance       = 0.555  (se = 0.0325)
#>   Likelihood ratio  = 5.07  on 2 df,  p = 0.07907
#>   Wald test         = 5.17  on 2 df,  p = 0.07534
#>   Score (logrank)   = 5.22  on 2 df,  p = 0.07358
#> =================================================================
#> 
```
