# Getting Started with BJM

## Overview

BJM fits a *backward joint model* of multivariate longitudinal outcomes
and time-to-event data, and uses it to make **dynamic predictions**:
given a patient’s longitudinal history up to some `prediction_time`,
predict their risk of an event (and, if desired, their future biomarker
values) conditional on survival to that point.

Fitting and prediction is a four-step pipeline:

1.  [`survivalSub()`](https://liwh0904.github.io/BJM/reference/survivalSub.md)
    fits the marginal survival sub-model (and, optionally, a
    competing-risks sub-model).
2.  [`longitudinalSub()`](https://liwh0904.github.io/BJM/reference/longitudinalSub.md)
    fits the longitudinal sub-model(s), one per biomarker.
3.  [`predictRisk()`](https://liwh0904.github.io/BJM/reference/predictRisk.md)
    combines the two fits to predict a patient’s event risk over a
    future horizon.
4.  [`predictLongitudinal()`](https://liwh0904.github.io/BJM/reference/predictLongitudinal.md)
    additionally predicts a future value of one or more specific
    biomarkers.

This vignette walks through all four steps on the `pbc3` data set that
ships with the package, then covers two practical questions that come up
once the mechanics are working: how to choose the `bandcount1`/
`bandcount2`/`bandcount3` numerical-integration tuning parameters, and
what happens when the pipeline functions are called with malformed
inputs.

``` r

library(BJM)
#> Loading required package: survival
data(pbc3)
```

`pbc3` is in long format: one row per longitudinal measurement, with
baseline covariates (`age`, `sex`, …) repeated on every row for a given
patient, and event-time variables (`years`, `status3`, …) that are also
constant within patient.

## Step 1: Fit the survival sub-model with `survivalSub()`

[`survivalSub()`](https://liwh0904.github.io/BJM/reference/survivalSub.md)
needs one row per patient, so we first drop the repeated longitudinal
rows. `form_marginal_surv` is a standard
[`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) formula
for the overall event. `form_conditional_cr` is optional: supply it when
there is more than one competing event type and you want to predict
event-specific risk (here, `status4` distinguishes the two causes among
patients who had *some* event).

``` r

data_survival_fitting <- pbc3[!duplicated(pbc3$id), ]

survival_fit_all <- survivalSub(
  data_survival_fitting,
  form_marginal_surv = Surv(years, status3) ~ age + sex,
  form_conditional_cr = status4 ~ years + age + sex
)

survival_fit_all
#> 
#> Call:
#> survivalSub(form_marginal_surv = Surv(years, status3) ~ age + sex,
#>             form_conditional_cr = status4 ~ years + age + sex)
#> 
#> Data Descriptives:
#>   Number of subjects        : 312
#>   Number of events          : 169
#>   Cause-1 events (CR model) : 29
#>   Cause-2 events (CR model) : 140
#> 
#> =================================================================
#>  Marginal Survival Sub-model  [Cox PH]
#> -----------------------------------------------------------------
#>  Formula: Surv(years, status3) ~ age + sex
#> 
#>          Coef exp(Coef)        SE      z p-value   
#> age  0.020411  1.020621  0.007584  2.691 0.00712 **
#> sex -0.497707  0.607923  0.207078 -2.403 0.01624 * 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#>   n = 312,  events = 169
#>   Concordance       = 0.570  (se = 0.0226)
#>   Likelihood ratio  = 14.20  on 2 df,  p = 0.0008263
#>   Wald test         = 15.29  on 2 df,  p = 0.0004776
#>   Score (logrank)   = 15.50  on 2 df,  p = 0.0004299
#> 
#> =================================================================
#>  Conditional Competing-Risks Sub-model  [Logistic GLM]
#> -----------------------------------------------------------------
#>  Formula: status4 ~ years + age + sex
#> 
#>                 Coef       SE      z  p-value    
#> (Intercept)  5.65622  1.70478  3.318 0.000907 ***
#> years       -0.02898  0.08919 -0.325 0.745242    
#> age         -0.15305  0.03076 -4.976  6.5e-07 ***
#> sex          0.02109  0.77150  0.027 0.978192    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#>   Null deviance     : 154.94  on 168 df
#>   Residual deviance : 115.05  on 165 df
#>   AIC : 123.05
#> =================================================================
```

## Step 2: Fit the longitudinal sub-models with `longitudinalSub()`

Each biomarker gets its own fixed-effects formula and random-effects
formula, supplied as same-length, same-order lists. `data_fit_all` is a
matching list of fitting data sets (one per biomarker) — or a single
`data.frame` that is reused for every biomarker, if all biomarkers are
fit on the same data. Here we fit two biomarkers, `serBilir` and
`albumin`, on the complete-case subset.

``` r

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

long_fit_all
#> 
#> Call:
#> longitudinalSub(M = 2 longitudinal outcomes)
#> 
#> =================================================================
#>  Longitudinal Sub-model(s)
#> -----------------------------------------------------------------
#>  [1] serBilir (continuous)
#>  Formula: serBilir ~ year + age + sex + (years) + (years) * year
#> 
#>                 Value        SE      t  p-value    
#> (Intercept)  2.424554  0.384604  6.304 5.13e-10 ***
#> year         0.442242  0.042957 10.295  < 2e-16 ***
#> age         -0.011675  0.005785 -2.018   0.0439 *  
#> sex          0.078449  0.167140  0.469   0.6390    
#> years       -0.204130  0.021557 -9.469  < 2e-16 ***
#> year:years  -0.027852  0.006112 -4.557 6.12e-06 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#>   Residual std. error (sigma): 0.4186
#>   n (subjects) = 169,  N (observations) = 872
#>   Log-likelihood: -772.87
#>   AIC: 1565.74,  BIC: 1613.45
#> 
#> -----------------------------------------------------------------
#>  [2] albumin (continuous)
#>  Formula: albumin ~ year + age + sex + (years) + (years) * year
#> 
#>                 Value        SE      t  p-value    
#> (Intercept)  3.578828  0.150261 23.817  < 2e-16 ***
#> year        -0.219541  0.026452 -8.300 5.36e-16 ***
#> age         -0.006625  0.002261 -2.930  0.00350 ** 
#> sex         -0.170236  0.064981 -2.620  0.00899 ** 
#> years        0.077939  0.008732  8.926  < 2e-16 ***
#> year:years   0.010955  0.003585  3.056  0.00233 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#>   Residual std. error (sigma): 0.3121
#>   n (subjects) = 169,  N (observations) = 872
#>   Log-likelihood: -386.29
#>   AIC: 792.58,  BIC: 840.29
#> 
#> =================================================================
#>  Multivariate Random-Effects Covariance Matrix (D)
#> -----------------------------------------------------------------
#>               (Intercept)_1  year_1 (Intercept)_2  year_2
#> (Intercept)_1        0.6012 -0.0231       -0.0473 -0.0070
#> year_1              -0.0231  0.0193        0.0121 -0.0050
#> (Intercept)_2       -0.0473  0.0121        0.0722 -0.0045
#> year_2              -0.0070 -0.0050       -0.0045  0.0050
#> =================================================================
```

## Optional: imputing interrupted follow-up before fitting

[`longitudinalSub()`](https://liwh0904.github.io/BJM/reference/longitudinalSub.md),
as used above, does complete-case analysis: for each biomarker it drops
rows missing that biomarker’s own covariates, and then keeps only the
subjects who have at least one non-missing observation of *every*
biomarker. With interrupted/irregular follow-up this throws away real
information – a subject who has `serBilir` at every visit but never had
`albumin` drawn is dropped from *both* models, not just the `albumin`
one.

[`imputeLongitudinal()`](https://liwh0904.github.io/BJM/reference/imputeLongitudinal.md)
is an optional preprocessing step that fills these gaps before
[`longitudinalSub()`](https://liwh0904.github.io/BJM/reference/longitudinalSub.md)
runs, assuming the missingness is at random (MAR) given the observed
covariates and biomarkers. Two backends are available, selected with
`method`:

- `method = "miwae"` (default) fits a MIWAE (Missing data
  Importance-Weighted AutoEncoder; Mattei & Frellsen, 2019) jointly
  across the supplied biomarkers and draws plausible completions via
  self-normalized importance resampling.
- `method = "diffusion"` fits a conditional denoising diffusion
  probabilistic model (DDPM; Ho, Jain & Abbeel, 2020) instead, trained
  with a masked score-matching loss, and imputes with a RePaint-style
  (Lugmayr et al., 2022) reverse-diffusion sampler that reproduces every
  observed value exactly and only draws the missing ones.

Both share the same call shape and return shape, so switching between
them is a one-argument change; the completed data can be passed straight
to
[`longitudinalSub()`](https://liwh0904.github.io/BJM/reference/longitudinalSub.md)
in place of the original `data_fit_all` either way. It requires the
optional `torch` package (listed in `Suggests`, not `Imports`), so the
rest of this section only runs if `torch` is installed.

``` r

torch_available <- requireNamespace("torch", quietly = TRUE) &&
  isTRUE(tryCatch(torch::torch_is_installed(), error = function(e) FALSE))
torch_available
#> [1] FALSE
```

To see it in action, artificially delete 12% of the `serBilir` and
`albumin` measurements in the fitting data at random – satisfying MAR by
construction, and mimicking a subject whose visit happened but one
biomarker wasn’t drawn that day:

``` r

set.seed(1)
data_fit_masked <- data_fit_all[[1]]
n <- nrow(data_fit_masked)
data_fit_masked$serBilir[sample.int(n, floor(0.12 * n))] <- NA
data_fit_masked$albumin[sample.int(n, floor(0.12 * n))] <- NA
```

``` r

imputed <- imputeLongitudinal(
  data_fit_masked, long_sub_fixed, long_sub_random,
  time_variable = "year", epochs = 80, seed = 1
)
```

`diagnostics` reports, per biomarker, the observed mean/sd alongside the
deep-generative and classical-baseline imputed mean/sd – useful for
judging whether the deep model’s completions look reasonable before
trusting them on your own data, since a deep generative model needs
enough data to fit reliably and is not automatically the better choice
at every sample size:

``` r

imputed$diagnostics$serBilir
```

Feeding the completed data into
[`longitudinalSub()`](https://liwh0904.github.io/BJM/reference/longitudinalSub.md)
retains at least as many subjects as running it on the masked data
directly, because a subject who was excluded only for missing one
biomarker at one visit is no longer excluded once that cell is filled
in:

``` r

count_retained <- function(data_fit_all, long_sub_fixed, id = "id") {
  if (is.data.frame(data_fit_all)) {
    data_fit_all <- rep(list(data_fit_all), length(long_sub_fixed))
  }
  ids_by_marker <- lapply(seq_along(long_sub_fixed), function(m) {
    d <- data_fit_all[[m]]
    complete <- d[rowSums(is.na(d[all.vars(long_sub_fixed[[m]])])) == 0, ]
    unique(complete[[id]])
  })
  Reduce(intersect, ids_by_marker)
}

length(count_retained(data_fit_masked, long_sub_fixed))
length(count_retained(imputed$data_fit_all, long_sub_fixed))
```

``` r

long_fit_imputed <- longitudinalSub(imputed$data_fit_all, long_sub_fixed, long_sub_random)
```

### A diffusion-based alternative

The same masked data can be completed with `method = "diffusion"`
instead – same arguments, a different generative model underneath:

``` r

imputed_diffusion <- imputeLongitudinal(
  data_fit_masked, long_sub_fixed, long_sub_random,
  time_variable = "year", method = "diffusion",
  epochs = 60, diffusion_steps = 40, seed = 1
)
imputed_diffusion$diagnostics$serBilir
```

``` r

long_fit_imputed_diffusion <- longitudinalSub(imputed_diffusion$data_fit_all,
                                               long_sub_fixed, long_sub_random)
length(count_retained(imputed_diffusion$data_fit_all, long_sub_fixed))
```

Neither backend is automatically the better choice: compare the
`diagnostics` above against the observed mean/sd, and, if in doubt,
compare downstream fixed-effect estimates from `long_fit_imputed` and
`long_fit_imputed_diffusion` against `long_fit_all` (the original
complete-case fit) before relying on either for inference.

### Propagating imputation uncertainty with `poolLongitudinalSub()`

`imputed` above used the default `impute = "single"`: `data_fit_all` is
filled with the *mean* of the completions drawn from the generative
model, giving one completed dataset. That is convenient, but it
understates uncertainty – fitting
[`longitudinalSub()`](https://liwh0904.github.io/BJM/reference/longitudinalSub.md)
on a single completed dataset treats the imputed values as if they were
observed data, so the resulting standard errors do not reflect not
actually knowing them.

`impute = "multiple"` instead returns `n_imputations` separately
completed datasets in `data_fit_all_list`, without averaging them
together. Fitting
[`longitudinalSub()`](https://liwh0904.github.io/BJM/reference/longitudinalSub.md)
once per completion and combining the fixed-effect estimates with
Rubin’s rules (Rubin, 1987) – implemented here by
[`poolLongitudinalSub()`](https://liwh0904.github.io/BJM/reference/poolLongitudinalSub.md)
– propagates the disagreement between completions into the final
standard errors, instead of discarding it:

``` r

imputed_multiple <- imputeLongitudinal(
  data_fit_masked, long_sub_fixed, long_sub_random,
  time_variable = "year", n_imputations = 5, impute = "multiple",
  epochs = 80, seed = 1
)
length(imputed_multiple$data_fit_all_list)
```

``` r

long_fit_all_list <- lapply(imputed_multiple$data_fit_all_list, function(d) {
  longitudinalSub(d, long_sub_fixed, long_sub_random)
})
pooled <- poolLongitudinalSub(long_fit_all_list)
pooled
```

The printed table’s `FMI` (fraction of missing information) column
reports, per coefficient, how much precision the missingness cost
relative to having observed everything – a coefficient with a high FMI
is one where the choice to impute (and how) matters more for that
specific estimate, even if the biomarker’s overall missingness rate
looks modest. `pooled$pooled` is a named list (one `data.frame` per
biomarker) with the pooled `estimate`, `std_error`, `df`, `statistic`,
`p_value`, `riv`, and `fmi` for programmatic use.

## Step 3: Predict event risk with `predictRisk()`

To predict for a specific patient, build a `data_predict_all` list (one
`data.frame` per biomarker, mirroring `data_fit_all`) containing only
that patient’s measurements *up to* `prediction_time` — later
measurements would not be available yet in a real prediction setting.
Below we predict, for patient 2, the risk of each event type within one
year of `year = 5`.

`survival_variable_all`/`survival_trans_function` describe how the raw
event-time variable is transformed for the integration grid; see
[`?predictRisk`](https://liwh0904.github.io/BJM/reference/predictRisk.md)
for details.

``` r

survival_variable_all <- list("Tyears1", "Tyears2", "Tyears3", "Tyears4")
survival_trans_function <- list(
  fun1 = function(x) abs(x - 1),
  fun2 = function(x) abs(x - 3),
  fun3 = function(x) abs(x - 5),
  fun4 = function(x) abs(x - 7)
)

data_raw_predict <- pbc3[pbc3$id == 2, ]
data_predict_all <- list(data_raw_predict, data_raw_predict)

risk <- predictRisk(
  data_predict_all, long_fit_all, survival_fit_all,
  prediction_time = 5, horizon = 1, time_variable = "year",
  survival_variable_all, survival_trans_function,
  bandcount1 = 10, bandcount2 = 20
)

risk
#> 
#> =================================================================
#>  Dynamic Prediction - Event Risk
#> -----------------------------------------------------------------
#>   Competing risks   : Yes
#>   Subjects          : 1
#> -----------------------------------------------------------------
#> 
#>  Subject Cause 1 Risk Cause 2 Risk Total Risk
#>       S1       0.0353       0.0016     0.0369
#> 
#> =================================================================
```

`survival_variable_all`/`survival_trans_function` almost always follow
the same convention shown above: variables named `"Tyears1"`,
`"Tyears2"`, …, each the absolute distance from a fixed cut point.
[`survivalTrans()`](https://liwh0904.github.io/BJM/reference/survivalTrans.md)
builds exactly that pair from a plain vector of cut points, so you do
not have to hand-write two matching parallel lists:

``` r

trans <- survivalTrans(c(1, 3, 5, 7))
identical(trans$survival_variable_all, survival_variable_all)
#> [1] TRUE
trans$survival_trans_function[[1]](2)
#> [1] 1
```

`risk_prob_1` and `risk_prob_2` are the predicted probabilities of
experiencing each of the two competing event types within the one-year
horizon, conditional on the patient’s longitudinal history and survival
to `prediction_time`. (`risk_prob_2` is `NULL` whenever
`survival_fit_all` was fit without `form_conditional_cr`.)

## Step 4: Predict a future biomarker value with `predictLongitudinal()`

[`predictLongitudinal()`](https://liwh0904.github.io/BJM/reference/predictLongitudinal.md)
answers a different question: not *whether* an event happens, but what a
specific biomarker’s value is likely to be at
`prediction_time + horizon`, conditional on survival. `bio_i` selects
the biomarker(s) by position in `long_fit_all` (`1` = `serBilir` here);
naming exactly one biomarker returns a single result, while naming more
than one (or leaving `bio_i = NULL` for every biomarker) returns a named
list of results instead, computing the shared, bio_i-independent part of
the pipeline only once and reusing it across every biomarker.

``` r

bio_pred <- predictLongitudinal(
  bio_i = 1, data_predict_all, long_fit_all, survival_fit_all,
  prediction_time = 5, horizon = 1, time_variable = "year",
  survival_variable_all, survival_trans_function,
  bandcount2 = 20, bandcount3 = 50
)

bio_pred$Y_predict
#> [1] 1.03628
```

`Y_predict` is the MAP (most likely) predicted value;
`Y_density`/`Y_all` give the full predicted density over a grid of
candidate values, which is what
[`predictPlot()`](https://liwh0904.github.io/BJM/reference/predictPlot.md)
visualizes (see
[`?predictPlot`](https://liwh0904.github.io/BJM/reference/predictPlot.md)).

## Choosing `bandcount1`, `bandcount2`, `bandcount3`

[`predictRisk()`](https://liwh0904.github.io/BJM/reference/predictRisk.md)
and
[`predictLongitudinal()`](https://liwh0904.github.io/BJM/reference/predictLongitudinal.md)
estimate their outputs by numerical integration over patient-specific
time and biomarker grids. The `bandcount*` arguments control how fine
those grids are:

- **`bandcount1`**
  ([`predictRisk()`](https://liwh0904.github.io/BJM/reference/predictRisk.md)
  only) is the number of grid points spanning the prediction window
  itself, from `prediction_time` to `prediction_time + horizon`. This is
  the numerator of the risk probability.
- **`bandcount2`** is the number of grid points spanning
  `[prediction_time, upper_bound]`, where `upper_bound` is set
  internally to twice the longest observed survival/censoring time among
  at-risk patients. This approximates “integrating out to infinity” for
  the denominator that normalizes the risk probability. A wider
  follow-up range needs a larger `bandcount2` to keep the grid spacing
  comparable.
- **`bandcount3`**
  ([`predictLongitudinal()`](https://liwh0904.github.io/BJM/reference/predictLongitudinal.md)
  only) is the number of points in the candidate-biomarker-value grid
  (`Y_all`) used to build the predicted density and locate its mode
  (`Y_predict`). It controls the resolution of the density curve, not a
  time integral.

Larger values give more accurate, smoother results at the cost of more
computation. There is no universal “correct” value because it depends on
how quickly the underlying hazard and biomarker trajectories change and
on the length of follow-up in your data — so the practical approach is a
**convergence check**: run the prediction once with the defaults, once
with every `bandcount*` doubled, and confirm the results barely move.

``` r

risk_default <- predictRisk(
  data_predict_all, long_fit_all, survival_fit_all,
  prediction_time = 5, horizon = 1, time_variable = "year",
  survival_variable_all, survival_trans_function,
  bandcount1 = 10, bandcount2 = 20
)

risk_doubled <- predictRisk(
  data_predict_all, long_fit_all, survival_fit_all,
  prediction_time = 5, horizon = 1, time_variable = "year",
  survival_variable_all, survival_trans_function,
  bandcount1 = 20, bandcount2 = 40
)

abs(risk_default$risk_prob_1 - risk_doubled$risk_prob_1)
#> [1] 0.003635748
abs(risk_default$risk_prob_2 - risk_doubled$risk_prob_2)
#> [1] 0.0001577366
```

If doubling the `bandcount*` values changes the result by more than you
can tolerate, keep doubling until it doesn’t; if it barely changes
anything (as above), the smaller, cheaper value is fine to use. The same
check applies to `bandcount3` for
[`predictLongitudinal()`](https://liwh0904.github.io/BJM/reference/predictLongitudinal.md).

## Friendly error messages

The pipeline functions validate their arguments before doing any model
fitting or numerical integration, so common mistakes fail fast with a
message that names the offending argument — instead of a cryptic error
from deep inside
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) or list
indexing. For example, passing a single `data.frame` instead of a list
to
[`predictRisk()`](https://liwh0904.github.io/BJM/reference/predictRisk.md):

``` r

predictRisk(
  data_predict_all[[1]], long_fit_all, survival_fit_all,
  prediction_time = 5, horizon = 1, time_variable = "year",
  survival_variable_all, survival_trans_function,
  bandcount1 = 10, bandcount2 = 20
)
#> 
#> =================================================================
#>  Dynamic Prediction - Event Risk
#> -----------------------------------------------------------------
#>   Competing risks   : Yes
#>   Subjects          : 1
#> -----------------------------------------------------------------
#> 
#>  Subject Cause 1 Risk Cause 2 Risk Total Risk
#>       S1       0.0353       0.0016     0.0369
#> 
#> =================================================================
```

or a formula referencing a column that doesn’t exist:

``` r

longitudinalSub(pbc3, serBilir ~ year + not_a_column, ~ year | id)
#> Error:
#> ! Variable(s) used in `long_sub_fixed[[1]]/long_sub_random[[1]]` not found in `data_fit_all[[1]]`: not_a_column.
```

## Where to go next

- [`?predictPlot`](https://liwh0904.github.io/BJM/reference/predictPlot.md)
  and
  [`?riskPlot`](https://liwh0904.github.io/BJM/reference/riskPlot.md)
  visualize the outputs of
  [`predictLongitudinal()`](https://liwh0904.github.io/BJM/reference/predictLongitudinal.md)/[`predictRisk()`](https://liwh0904.github.io/BJM/reference/predictRisk.md)
  for a single patient.
- [`?cmtPlot`](https://liwh0904.github.io/BJM/reference/cmtPlot.md)
  plots observed longitudinal trajectories stratified by eventual
  outcome, useful for checking whether a biomarker looks informative
  before fitting.
- [`?imputeLongitudinal`](https://liwh0904.github.io/BJM/reference/imputeLongitudinal.md)
  fills in missing longitudinal biomarker values due to interrupted
  follow-up before
  [`longitudinalSub()`](https://liwh0904.github.io/BJM/reference/longitudinalSub.md)
  runs, instead of relying on complete-case analysis (see “Optional:
  imputing interrupted follow-up before fitting” above).
- [`?poolLongitudinalSub`](https://liwh0904.github.io/BJM/reference/poolLongitudinalSub.md)
  combines
  [`longitudinalSub()`](https://liwh0904.github.io/BJM/reference/longitudinalSub.md)
  fits across the multiple completed datasets from
  `imputeLongitudinal(..., impute = "multiple")` with Rubin’s rules, so
  imputation uncertainty is reflected in the reported standard errors
  (see “Propagating imputation uncertainty with
  [`poolLongitudinalSub()`](https://liwh0904.github.io/BJM/reference/poolLongitudinalSub.md)”
  above).
