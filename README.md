# BJM

<!-- badges: start -->
<!-- badges: end -->

BJM fits a *backward joint model* of multivariate longitudinal outcomes
and time-to-event data, and uses it to make **dynamic predictions**:
given a patient's longitudinal history up to some `prediction_time`,
predict their risk of an event (and, if desired, their future biomarker
values) conditional on survival to that point. The backward joint model
formulation keeps fitting and prediction fast, even with large sample
sizes and many longitudinal variables, and accommodates irregularly
measured longitudinal data and competing risks.

## Installation

```r
install.packages("BJM")
```

To install the development version from source (e.g. a local clone or
tarball):

```r
# install.packages("devtools")
devtools::install_local("path/to/BJM")
```

## Quick start

Fitting and prediction is a four-step pipeline: fit a survival
sub-model, fit a longitudinal sub-model per biomarker, then combine
them to predict event risk and/or future biomarker values.

```r
library(BJM)
data(pbc3)

## 1. Fit the survival sub-model (one row per patient)
data_survival_fitting <- pbc3[!duplicated(pbc3$id), ]

survival_fit_all <- survivalSub(
  data_survival_fitting,
  form_marginal_surv = Surv(years, status3) ~ age + sex,
  form_conditional_cr = status4 ~ years + age + sex
)

## 2. Fit the longitudinal sub-model(s), one per biomarker
long_sub_fixed <- list(
  "serBilir" = serBilir ~ year + age + sex + (years) + (years) * year,
  "albumin"  = albumin  ~ year + age + sex + (years) + (years) * year
)
long_sub_random <- list(
  "serBilir" = ~ year | id,
  "albumin"  = ~ year | id
)
data_fit_all <- list(pbc3[pbc3$status3 == 1, ], pbc3[pbc3$status3 == 1, ])

long_fit_all <- longitudinalSub(data_fit_all, long_sub_fixed, long_sub_random)

## 3. Predict event risk for a patient, using only their history up to
##    prediction_time
survival_variable_all <- list("Tyears1", "Tyears2", "Tyears3", "Tyears4")
survival_trans_function <- list(
  fun1 = function(x) abs(x - 1),
  fun2 = function(x) abs(x - 3),
  fun3 = function(x) abs(x - 5),
  fun4 = function(x) abs(x - 7)
)

data_raw_predict <- pbc3[pbc3$id == 2, ]
data_predict_all <- list(data_raw_predict, data_raw_predict)

risk <- dynamicPrediction(
  data_predict_all, long_fit_all, survival_fit_all,
  prediction_time = 5, horizon = 1, time_variable = "year",
  survival_variable_all, survival_trans_function,
  bandcount1 = 10, bandcount2 = 20
)
risk

## 4. Predict a future biomarker value conditional on survival
bio_pred <- dynamicPredictionBio(
  bio_i = 1, data_predict_all, long_fit_all, survival_fit_all,
  prediction_time = 5, horizon = 1, time_variable = "year",
  survival_variable_all, survival_trans_function,
  bandcount2 = 20, bandcount3 = 50
)
bio_pred$Y_predict
```

`survival_variable_all`/`survival_trans_function` above follow a common
convention (`"Tyears1"`, `"Tyears2"`, ..., each the absolute distance from a
fixed cut point); `survivalTrans(c(1, 3, 5, 7))` builds that same pair for
you instead of hand-writing two matching lists. And every `data_*_all`
argument shown above (`data_fit_all`, `data_predict_all`, ...) also accepts
a single bare `data.frame` — reused for every biomarker — instead of a
repeated list, when all biomarkers share the same measurement data.

`predictPlot()` and `riskPlot()` visualize these predictions for a
single patient; `cmtPlot()` plots observed longitudinal trajectories
stratified by eventual outcome.

For the full walkthrough — including how to choose the `bandcount1`/
`bandcount2`/`bandcount3` numerical-integration tuning parameters via a
convergence check, and what the input-validation error messages look
like — see `vignette("BJM-intro", package = "BJM")`.
