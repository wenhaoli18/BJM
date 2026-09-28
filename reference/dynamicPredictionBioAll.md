# Predict multiple future biomarker values from fitted sub-models

**Internal multi-biomarker engine** behind
[`predictLongitudinal`](https://wenhaoli18.github.io/BJM/reference/predictLongitudinal.md)
– call
[`predictLongitudinal()`](https://wenhaoli18.github.io/BJM/reference/predictLongitudinal.md)
directly instead (it dispatches here automatically when `bio_i` names
more than one biomarker, or is left `NULL`).

Batch counterpart to
[`dynamicPredictionBio`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md):
predicts every requested biomarker in `bio_i` from the same fitted
longitudinal/survival sub-models and the same `data_predict_all`, but
computes the bio_i-**independent** part of the pipeline (restricting to
at-risk patients, the survival-side integration grid, and the
denominator conditional density – see
[`compute_bio_shared_step()`](https://wenhaoli18.github.io/BJM/reference/compute_bio_shared_step.md))
only **once** and reuses it across every biomarker, instead of
recomputing it once per biomarker the way calling
[`dynamicPredictionBio()`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md)
in a loop would. This matters most under the Gaussian-copula path (see
[`longitudinalSubCopula`](https://wenhaoli18.github.io/BJM/reference/longitudinalSubCopula.md)),
where that denominator involves
[`mvtnorm::pmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/pmvnorm.html)
Monte-Carlo evaluations that are otherwise the dominant cost of a
multi-biomarker prediction run.

`bandcount2` (controlling the shared survival-integration grid) and
`bandcount3` (controlling each biomarker's own candidate-value grid) are
auto-tuned separately when left at their default `"auto"`: if
`bandcount2 = "auto"`, it is resolved **once**, using a single
representative biomarker (the first one in `bio_i`), via the same
doubling-until-stable check
[`dynamicPredictionBio`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md)
uses; the shared step is then built once at that resolved value. If
`bandcount3 = "auto"`, it is then resolved independently for **every**
biomarker in `bio_i` (their candidate-value grids need not converge at
the same resolution), reusing the once-computed shared step for every
doubling round rather than rebuilding it – except for any **ordinal**
biomarker, for which `bandcount3` tuning is always skipped (its
candidate grid is fixed at its category count; see
[`dynamicPredictionBio`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md))
and `"bandcount3"` is recorded as `NA` for that biomarker. See Details
in
[`dynamicPredictionBio`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md)
for the general auto-bandcount rationale.

## Usage

``` r
dynamicPredictionBioAll(
  bio_i = NULL,
  data_predict_all,
  long_fit_all,
  survival_fit_all,
  prediction_time,
  horizon,
  time_variable,
  survival_variable_all,
  survival_trans_function,
  bandcount2 = "auto",
  bandcount3 = "auto"
)
```

## Arguments

- bio_i:

  Integer vector of biomarkers to predict. May include continuous and/or
  ordinal biomarkers (see
  [`dynamicPredictionBio`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md)).
  Defaults to `NULL`, meaning every biomarker in `long_fit_all`.

- data_predict_all:

  See
  [`dynamicPredictionBio`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md).

- long_fit_all:

  See
  [`dynamicPredictionBio`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md).

- survival_fit_all:

  See
  [`dynamicPredictionBio`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md).

- prediction_time:

  See
  [`dynamicPredictionBio`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md).

- horizon:

  See
  [`dynamicPredictionBio`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md).

- time_variable:

  See
  [`dynamicPredictionBio`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md).

- survival_variable_all:

  See
  [`dynamicPredictionBio`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md).

- survival_trans_function:

  See
  [`dynamicPredictionBio`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md).

- bandcount2:

  See
  [`dynamicPredictionBio`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md).

- bandcount3:

  See
  [`dynamicPredictionBio`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md).

## Value

A named list of `"dynamicPredictionBio.BJM"` objects (see
[`dynamicPredictionBio`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md)),
one per requested biomarker, named by that biomarker's response-variable
name; with attributes `"bandcount2"` (the single resolved/used
`bandcount2`) and `"bandcount3"` (a named numeric vector of the
resolved/used `bandcount3` for each biomarker, or `NA` for an ordinal
biomarker). Classed `"dynamicPredictionBioAll.BJM"`.

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

Y_predict_all = dynamicPredictionBioAll(bio_i = NULL, data_predict_all, long_fit_all,
                                        survival_fit_all, prediction_time = 3,
                                        horizon = 3, time_variable = "year",
                                        survival_variable_all, survival_trans_function,
                                        bandcount2 = 40, bandcount3 = 400)
#> Error in dynamicPredictionBioAll(bio_i = NULL, data_predict_all, long_fit_all,     survival_fit_all, prediction_time = 3, horizon = 3, time_variable = "year",     survival_variable_all, survival_trans_function, bandcount2 = 40,     bandcount3 = 400): could not find function "dynamicPredictionBioAll"

# }
```
