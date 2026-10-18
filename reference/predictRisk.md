# Dynamic prediction function for future event risk

Combines a fitted longitudinal sub-model
([`longitudinalSub`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md))
and survival sub-model
([`survivalSub`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md))
into a backward joint model prediction: for each at-risk subject at
`prediction_time`, the subject's observed longitudinal history up to
that time is used to update their individual random effects (empirical
Bayes), which are then integrated against the survival sub-model's
hazard to give the predicted probability of experiencing the event
within `(prediction_time, prediction_time + horizon]`, conditional on
being event-free at `prediction_time`. The numerator (event within the
window) is evaluated on a grid over the window whose resolution is
controlled by `bandcount1`; the denominator (survival to
`prediction_time`, used for normalization) is that same window integral
plus the integral beyond the window, evaluated on a grid controlled by
`bandcount2`; see Details.

The prediction is conditional on the longitudinal history observed up to
`prediction_time`: rows of `data_predict_all` whose `time_variable` is
later than `prediction_time` are dropped, with a warning, before
predicting.

A patient with no biomarker measurement up to `prediction_time` (rows
whose biomarker values are all `NA`, or only rows after
`prediction_time`, e.g. a `prediction_time` before the first visit) has
no history to condition on. The longitudinal factor \\f(Y \mid T)\\ of
the model is then 1, and the patient's risk is that of the survival
sub-model alone, given the covariates in their earliest row; a warning
names these patients. The other patients' predictions do not change. A
patient with measurements of some biomarkers but none of another is
still not predicted (`NA`, with a warning).

## Usage

``` r
predictRisk(
  data_predict_all,
  long_fit_all,
  survival_fit_all,
  prediction_time,
  horizon,
  time_variable,
  survival_variable_all,
  survival_trans_function,
  bandcount1 = "auto",
  bandcount2 = "auto",
  n_cores = 1
)
```

## Arguments

- data_predict_all:

  This involves a collection of `data.frame` objects for dynamic
  prediction, each corresponding to a distinct longitudinal outcome.
  These data frames should contain the variables specified in
  `long_sub_fixed` and `long_sub_random`. Utilizing a list structure
  allows for the incorporation of multiple longitudinal outcomes, each
  potentially following different measurement protocols. In instances
  where all longitudinal outcomes are recorded at identical time points
  across patients, a singular `data.frame` object may be used in a
  `list`. Alternatively, a single bare `data.frame` (not wrapped in a
  list) may be supplied directly; it is then reused for every
  longitudinal outcome. It is presumed that each data frame is
  structured in a long format.

- long_fit_all:

  Outputs from the model fitting process using the `nlme` package,
  encompassing the results and parameters obtained from the analysis.

- survival_fit_all:

  Results and parameters generated from the model fitting procedure,
  utilizing the `coxph` function. These outputs include the
  comprehensive findings and variables derived from the analysis.

- prediction_time:

  Time used to make the prediction.

- horizon:

  Prediction horizon: the length of the window after `prediction_time`,
  `>= 0`.

- time_variable:

  The name of time variable in linear mixed model.

- survival_variable_all:

  The name of the transformed time-to-event outcomes variable. The
  longitudinal sub-model must use the event time, directly or through
  one of these transformations.

- survival_trans_function:

  The transformation function used for time-to-event outcomes, in the
  order of `survival_variable_all`. They are checked against the columns
  the longitudinal sub-model was fit on, and a mismatch (e.g. two
  functions in the wrong order) is an error.

- bandcount1:

  The number of grid points spanning the prediction window, from
  `prediction_time` to `prediction_time + horizon` (the numerator of the
  risk probability). Larger values give a more accurate but slower
  estimate. Defaults to `"auto"` (see Details).

- bandcount2:

  The number of grid points spanning
  `[prediction_time + horizon, upper_bound]`, where `upper_bound` is set
  internally as the earliest time by which every at-risk patient's
  model-based probability of still being event-free (given event-free at
  `prediction_time + horizon`) has dropped below `1e-4`; this
  approximates integrating out to infinity for the part of the
  denominator beyond the prediction window. A wider follow-up range
  needs a larger `bandcount2` to keep the grid spacing comparable.
  Defaults to `"auto"` (see Details).

- n_cores:

  Number of CPU cores. With an ordinal biomarker the prediction's cost
  is per patient, and the patients are spread over `n_cores` cores with
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html);
  this is only parallel on Unix-like systems (Linux, macOS), and runs
  serially on Windows. The result is exactly the same for any `n_cores`.
  With continuous biomarkers only, predictions are fast and `n_cores`
  has no effect. Defaults to `1`.

## Value

An object of class `"predictRisk.BJM"`, a named list with elements:

- risk_prob_1:

  A vector of dynamically predicted probabilities, one per at-risk
  patient, of experiencing the (first) event within the prediction
  horizon, named by patient id. All `0` when `horizon = 0`.

- risk_prob_2:

  When `survival_fit_all` was fit with competing risks, a vector of
  dynamically predicted probabilities, one per at-risk patient, of
  experiencing the competing event within the prediction horizon, named
  by patient id (all `0` when `horizon = 0`). `NULL` when there is no
  competing risk.

- prob_undetected_1, prob_undetected_2:

  Only for an interval-censored `survival_fit_all` (see below): the
  probability that the (first / competing) event already happened
  between the patient's last visit and `prediction_time`, undetected.
  `prob_undetected_2` is `NULL` without competing risks.

**Interval-censored fits.** An event is then only detected at a visit,
so at `prediction_time` \\s\\ a patient is known to be event-free only
up to their last visit \\V\\ (the last time in their history up to
\\s\\). `risk_prob_1` is then \\P(s \< T \le s + horizon \mid T \> V,
history)\\ and `prob_undetected_1` is \\P(V \< T \le s \mid T \> V,
history)\\ (with competing risks, each `_1`/`_2` value is that
probability jointly with the event type); with \\V = s\\ they reduce to
the usual conditioning on \\T \> s\\. Only patients still at risk at
`prediction_time` are predicted: a patient whose recorded survival time
is before `prediction_time` is left out (one with a missing survival
time is kept). The names show which patients each value belongs to. The
numerical integration grids are shared by all patients predicted in one
call, so a patient's value can differ slightly (within the grid's
discretization error, which shrinks as `bandcount1`/`bandcount2` grow)
depending on which other patients are predicted alongside it.

## Details

There is no universal correct value for `bandcount1`/`bandcount2`: as a
practical check, double both and confirm the resulting risk
probabilities barely change; if they do, keep doubling. By default
(`bandcount1 = "auto"`, `bandcount2 = "auto"`), this doubling check is
done for you: starting from small built-in values, both are doubled
together, and the result is compared to the previous round, until the
largest relative change in the risk probabilities drops below 1%, or 2
doublings have been tried (so at most 3 calls' worth of work). If it
still has not converged by then, a warning reports this and the result
at the largest value tried is returned anyway (not an error), so this
never silently loops for an unbounded amount of time. Pass an explicit
number for either argument to skip auto-tuning it and use a fixed value
instead (as in previous package versions), or call
[`checkBandcountConvergence()`](https://wenhaoli18.github.io/BJM/reference/checkBandcountConvergence.md)
directly for more control over the tolerance and doubling count. See
also
[`vignette("BJM-intro", package = "BJM")`](https://wenhaoli18.github.io/BJM/articles/BJM-intro.md)
for a worked example.

The denominator integrates over every event time after
`prediction_time`, including times beyond the last follow-up time in the
data
[`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md)
was fit on. There, the baseline cumulative hazard is extrapolated
linearly (a constant hazard), and the longitudinal sub-model's mean is
evaluated at event times it was never fit on. When much of an at-risk
patient's survival probability lies beyond the last follow-up time (e.g.
a prediction late in follow-up), the prediction depends on this
extrapolation; in `pbc3` it moved the risks checked by less than one
percentage point.

For a fit with an ordinal biomarker, a patient's ordinal measurements
contribute the probability of their observed categories. With a single
ordinal random effect (e.g. one ordinal biomarker with a random
intercept) it is computed by Gauss-Hermite quadrature, and is exact;
with more, by Monte Carlo
([`mvtnorm::pmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/pmvnorm.html)),
with an error around the fourth significant digit. Each patient's Monte
Carlo draws use a fixed seed derived from the patient's id, so the
results are the same on every run and for any `n_cores`, without
[`set.seed()`](https://rdrr.io/r/base/Random.html), and the caller's
random-number state is left unchanged.

## Examples

``` r

# \donttest{
data(pbc3)

data_survival_fitting =  pbc3[!duplicated(pbc3$id), ]

form_marginal_surv = Surv(years, status3) ~ age + sex
form_conditional_cr = NULL

survival_fit_all = survivalSub(data_survival_fitting, form_marginal_surv, 
                               form_conditional_cr)

long_sub_fixed = list(
  "long1" = serBilir ~ year + age + sex +  (years) + (years) * year,  
  "long2" = prothrombin ~ year + age + sex + (years) + (years) * year,  
  "long3" = albumin ~ year + age + age * year + sex + (years) + (years) * year,  
  "long4" = alkaline ~ year + age + sex + (years) + (years) * year, 
  "long5" = SGOT ~ year + age + sex + (years) + (years) * year, 
  "long6" = platelets ~ year + age + sex + (years)  + (years) * year)

long_sub_random =list(
  "long1" =  ~ year| id,   
  "long2" =  ~ year| id,    
  "long3" =  ~ year| id,    
  "long4" =  ~ year| id,    
  "long5" =  ~ year| id,    
  "long6" =  ~ year| id)

survival_variable_all = list(
  "Tyears1",  "Tyears2", "Tyears3", "Tyears4"
)

survival_trans_function = list(
  fun1 = function(x){abs(x - 1)}, 
  fun2 = function(x){abs(x - 3)}, 
  fun3 = function(x){abs(x - 5)}, 
  fun4 = function(x){abs(x - 7)}
)

# Complete case analysis
data_fit_all = list()
for(i in seq_len(length(long_sub_fixed))){
  data_fit_all[[i]] = pbc3[pbc3$status3 == 1, ]
}

# fitting longitudinal submodel
long_fit_all = longitudinalSub(data_fit_all, long_sub_fixed, long_sub_random)

i_PID = 2
data.raw.predict.1 = pbc3[pbc3$id == i_PID, ]

data_predict_all = list()
for(i in seq_len(length(long_sub_fixed))){
  data_predict_all[[i]] = data.raw.predict.1[data.raw.predict.1$year <= 3,]
}

# predict risk probability
risk.prob = predictRisk(data_predict_all, long_fit_all, survival_fit_all,
                              prediction_time = 3,
                              horizon = 3, time_variable = "year",
                              survival_variable_all, survival_trans_function,
                              bandcount1 = 10, bandcount2 = 10)

# poly() in its default orthogonal mode, splines::ns()/bs(), and factor()
# in a continuous biomarker's long_sub_fixed produce correct dynamic
# predictions (see ?longitudinalSub) -- including
# when a patient has only a single longitudinal observation to condition
# on -- because the basis/contrasts fit on the full training data are
# cached (via each biomarker's terms object and long_fit_all$xlevels) and
# reused here at prediction time, instead of being recomputed from that
# patient's small per-prediction-time slice:
long_sub_fixed_nonlinear = list(
  "long1" = serBilir ~ poly(year, 2) + age + sex + years,
  "long2" = albumin ~ year + age + sex + years)
long_sub_random_nonlinear = list("long1" = ~ year | id, "long2" = ~ year | id)
long_fit_nonlinear = longitudinalSub(pbc3[pbc3$status3 == 1, ],
                                     long_sub_fixed_nonlinear,
                                     long_sub_random_nonlinear)

data_predict_normal = data.raw.predict.1[data.raw.predict.1$year <= 3, ]
data_predict_sparse = data.raw.predict.1[1, ]

risk.prob.normal = predictRisk(data_predict_normal, long_fit_nonlinear,
                                     survival_fit_all, prediction_time = 3,
                                     horizon = 3, time_variable = "year",
                                     survival_variable_all, survival_trans_function,
                                     bandcount1 = 10, bandcount2 = 10)
risk.prob.sparse = predictRisk(data_predict_sparse, long_fit_nonlinear,
                                     survival_fit_all, prediction_time = 3,
                                     horizon = 3, time_variable = "year",
                                     survival_variable_all, survival_trans_function,
                                     bandcount1 = 10, bandcount2 = 10)
# both give a sane, non-degenerate risk_prob_1 (not 0, no error)
risk.prob.normal$risk_prob_1
#>        2 
#> 0.089629 
risk.prob.sparse$risk_prob_1
#>         2 
#> 0.1202103 

# }
```
