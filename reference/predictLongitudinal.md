# Dynamic prediction function for future longitudinal outcomes

Companion to
[`predictRisk`](https://liwh0904.github.io/BJM/reference/predictRisk.md),
using the same fitted longitudinal
([`longitudinalSub`](https://liwh0904.github.io/BJM/reference/longitudinalSub.md))
and survival
([`survivalSub`](https://liwh0904.github.io/BJM/reference/survivalSub.md))
sub-models, but instead of an event-risk probability this returns a
predictive density for a future value of one or more chosen biomarkers
(`bio_i`) at `prediction_time + horizon`, conditional on the subject's
observed longitudinal history up to `prediction_time` and on being
event-free at `prediction_time`.

This is the single, unified entry point for biomarker-value prediction:

- If `bio_i` names **exactly one** biomarker, this predicts just that
  biomarker and returns a single `"dynamicPredictionBio.BJM"` object.

- If `bio_i` names **more than one** biomarker, or is left at its
  default `NULL` (meaning every biomarker in `long_fit_all`), this
  instead predicts every requested biomarker and returns a named list of
  `"dynamicPredictionBio.BJM"` objects (classed
  `"dynamicPredictionBioAll.BJM"`), computing the bio_i-**independent**
  part of the pipeline (restricting to at-risk patients, the
  survival-side integration grid, and the denominator conditional
  density) only **once** and reusing it across every biomarker, instead
  of recomputing it once per biomarker. This matters most under the
  Gaussian-copula path (see
  [`longitudinalSubCopula`](https://liwh0904.github.io/BJM/reference/longitudinalSubCopula.md)),
  where that denominator involves
  [`mvtnorm::pmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/pmvnorm.html)
  Monte-Carlo evaluations that are otherwise the dominant cost of a
  multi-biomarker prediction run.

See Value for the exact return shape in each case.

`bio_i` may refer to either **continuous** or **ordinal** biomarkers
(see
[`longitudinalSubCopula`](https://liwh0904.github.io/BJM/reference/longitudinalSubCopula.md)),
and a single call may mix both. For an ordinal biomarker,
`Y_all`/`Y_predict` hold integer *category codes* (`1:K`, in the fitted
factor's [`levels()`](https://rdrr.io/r/base/levels.html)/threshold
order) rather than a numeric grid, with `Y_all` additionally carrying a
`"category_labels"` attribute giving the matching level-label strings;
`bandcount3` is ignored for that biomarker, since its candidate grid is
fixed at its category count (see Details).

The time values in the prediction data subset must be less than the
specified `prediction_time` which is the prediction time. The time
points for longitudinal repeated measurements must not surpass the
prediction time.

`bandcount2` (controlling the shared survival-integration grid) and
`bandcount3` (controlling each biomarker's own candidate-value grid) are
auto-tuned separately when left at their default `"auto"`, via the
doubling-until-stable check described in Details. When `bio_i` names
more than one biomarker, `bandcount2` is resolved **once** using a
single representative biomarker (the first one in `bio_i`), and
`bandcount3` is then resolved independently for **every** requested
biomarker (their candidate-value grids need not converge at the same
resolution) – except for any **ordinal** biomarker, for which
`bandcount3` tuning is always skipped (see above) and `"bandcount3"` is
recorded as `NA` for that biomarker.

## Usage

``` r
predictLongitudinal(
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
  ordinal biomarkers (see Details). Defaults to `NULL`, meaning every
  biomarker in `long_fit_all`. Naming exactly one biomarker returns a
  single result rather than a list of one – see Value.

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

  Time used to make the prediction

- horizon:

  Prediction horizon

- time_variable:

  The name of time variable in linear mixed model.

- survival_variable_all:

  The name of the transformed time-to-event outcomes variable.

- survival_trans_function:

  The transformation function used for time-to-event outcomes, in the
  order of `survival_variable_all`.

- bandcount2:

  The number of grid points spanning `[prediction_time, upper_bound]`,
  where `upper_bound` is set internally to twice the longest observed
  survival/censoring time among at-risk patients; this approximates
  integrating out to infinity for the denominator that normalizes the
  predicted density. A wider follow-up range needs a larger `bandcount2`
  to keep the grid spacing comparable. Defaults to `"auto"` (see
  Details).

- bandcount3:

  The number of points in the candidate-biomarker-value grid (`Y_all`)
  used to build the predicted density (`Y_density`) and locate its mode
  (`Y_predict`). This controls the resolution of the density curve, not
  a time integral; increase it if the density looks jagged or
  `Y_predict` jumps erratically between nearby grid points. Defaults to
  `"auto"` (see Details). Ignored for an ordinal biomarker (see Details)
  – that candidate grid is always the biomarker's fixed category count.

## Value

If `bio_i` names exactly one biomarker: an object of class
`"dynamicPredictionBio.BJM"`, a named list with elements:

- Y_predict:

  A vector, one entry per patient, giving the MAP (most likely)
  predicted value of the biomarker at `prediction_time + horizon`. For
  an ordinal biomarker, this is an integer category code (see Details),
  not a raw value.

- Y_density:

  A probability matrix whose rows correspond to the candidate biomarker
  values in `Y_all` and whose columns correspond to individual patients;
  each entry is the dynamically predicted density of the biomarker
  taking that value. For an ordinal biomarker, each row is instead that
  category's predicted probability.

- Y_all:

  The grid of candidate biomarker values used to build `Y_density`. For
  an ordinal biomarker, this is the vector of integer category codes
  `1:K`, carrying a `"category_labels"` attribute with the matching
  level-label strings (see Details).

Otherwise (`bio_i` names more than one biomarker, or is `NULL`): a named
list of such objects, one per requested biomarker, named by that
biomarker's response-variable name; with attributes `"bandcount2"` (the
single resolved/used `bandcount2`) and `"bandcount3"` (a named numeric
vector of the resolved/used `bandcount3` for each biomarker, or `NA` for
an ordinal biomarker). Classed `"dynamicPredictionBioAll.BJM"`.

## Details

There is no universal correct value for `bandcount2`/`bandcount3`: as a
practical check, double both and confirm the results barely change; if
they do, keep doubling. By default (`bandcount2 = "auto"`,
`bandcount3 = "auto"`), this doubling check is done for you: starting
from small built-in values, both are doubled together, and the result is
compared to the previous round, until the largest relative change in
`Y_predict` drops below 1%, or 2 doublings have been tried (so at most 3
calls' worth of work). If it still has not converged by then, a warning
reports this and the result at the largest value tried is returned
anyway (not an error), so this never silently loops for an unbounded
amount of time. Pass an explicit number for either argument to skip
auto-tuning it and use a fixed value instead, or call
[`checkBandcountConvergence()`](https://liwh0904.github.io/BJM/reference/checkBandcountConvergence.md)
directly for more control over the tolerance and doubling count. See
also
[`vignette("BJM-intro", package = "BJM")`](https://liwh0904.github.io/BJM/articles/BJM-intro.md)
for a worked example.

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

# A single biomarker -> a single dynamicPredictionBio.BJM result
Y_predict = predictLongitudinal(bio_i = 1, data_predict_all, long_fit_all,
                                survival_fit_all, prediction_time = 3,
                                horizon = 3, time_variable = "year",
                                survival_variable_all, survival_trans_function,
                                bandcount2 = 40, bandcount3 = 400)

# Multiple biomarkers (or bio_i = NULL for all of them) -> a named list
Y_predict_all = predictLongitudinal(bio_i = NULL, data_predict_all, long_fit_all,
                                    survival_fit_all, prediction_time = 3,
                                    horizon = 3, time_variable = "year",
                                    survival_variable_all, survival_trans_function,
                                    bandcount2 = 40, bandcount3 = 400)

# }
```
