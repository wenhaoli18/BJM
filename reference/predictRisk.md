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
being event-free at `prediction_time`. The integrals in both the
numerator (event within the window) and denominator (survival to
`prediction_time`, used for normalization) are evaluated on numerical
grids whose resolution is controlled by `bandcount1`/ `bandcount2`; see
Details.

The prediction is conditional on the longitudinal history observed up to
`prediction_time`: rows of `data_predict_all` whose `time_variable` is
later than `prediction_time` are dropped, with a warning, before
predicting.

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
  bandcount2 = "auto"
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

  The name of the transformed time-to-event outcomes variable.

- survival_trans_function:

  The transformation function used for time-to-event outcomes, in the
  order of `survival_variable_all`.

- bandcount1:

  The number of grid points spanning the prediction window, from
  `prediction_time` to `prediction_time + horizon` (the numerator of the
  risk probability). Larger values give a more accurate but slower
  estimate. Defaults to `"auto"` (see Details).

- bandcount2:

  The number of grid points spanning `[prediction_time, upper_bound]`,
  where `upper_bound` is set internally as the earliest time by which
  every at-risk patient's model-based probability of still being
  event-free (given event-free at `prediction_time`) has dropped below
  `1e-4`; this approximates integrating out to infinity for the
  denominator that normalizes the risk probability. A wider follow-up
  range needs a larger `bandcount2` to keep the grid spacing comparable.
  Defaults to `"auto"` (see Details).

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

Only patients still at risk at `prediction_time` are predicted: a
patient whose recorded survival time is before `prediction_time` is left
out (one with a missing survival time is kept). The names show which
patients each value belongs to. The numerical integration grids are
shared by all patients predicted in one call, so a patient's value can
differ slightly (within the grid's discretization error, which shrinks
as `bandcount1`/`bandcount2` grow) depending on which other patients are
predicted alongside it.

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

For a fit with an ordinal biomarker, each ordinal measurement
contributes a multivariate normal probability computed by Monte Carlo
([`mvtnorm::pmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/pmvnorm.html)),
so repeated calls differ slightly (around the fifth significant digit);
call [`set.seed()`](https://rdrr.io/r/base/Random.html) first for
exactly reproducible results.

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
#>          2 
#> 0.09646473 
risk.prob.sparse$risk_prob_1
#>        2 
#> 0.127205 

# }
```
