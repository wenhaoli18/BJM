# Fit a multivariate longitudinal sub-model

Fits one linear mixed-effects model per longitudinal biomarker
separately (via [`lme`](https://rdrr.io/pkg/nlme/man/lme.html)), keeping
each biomarker's own fixed-effects and residual-variance estimates from
that separate fit. The `M` separate fits are then combined into a single
multivariate model by re-estimating the full random-effects
variance-covariance matrix jointly across all `M` biomarkers, via an EM
algorithm initialized from the block-diagonal covariance implied by the
separate fits – so correlation between biomarkers' random effects is
captured, even though the fixed effects and residual variances
themselves are not re-estimated jointly and remain exactly what each
biomarker's own `lme()` fit produced. This lets you assess the effect of
predictor variables on several longitudinal outcomes at once while
accounting for both population-level (fixed) effects and subject-level
(random) variability, and for how the outcomes covary within a subject.
The result is one of the two sub-models – together with
[`survivalSub`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md)
– that
[`predictRisk`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)/[`dynamicPredictionBio`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md)
combine to produce dynamic risk/biomarker predictions.

## Usage

``` r
longitudinalSub(
  data_fit_all,
  long_sub_fixed,
  long_sub_random,
  biomarker_type = NULL
)
```

## Arguments

- data_fit_all:

  This process requires a set of `data.frame` objects designated for
  model fitting, with each `data.frame` representing a separate
  longitudinal outcome. These `data.frame` objects must include the
  variables identified in `long_sub_fixed` and `long_sub_random`. The
  use of a `list` arrangement facilitates the inclusion of various
  longitudinal outcomes, which may adhere to different measurement
  protocols. When all longitudinal outcomes are measured at the same
  time points for every patient, a single `data.frame` object can be in
  a list. It is assumed that every `data.frame` is organized in a long
  format.

- long_sub_fixed:

  This refers to a collection of formulas detailing the fixed effects
  portion for each longitudinal outcome. On the left side of each
  formula, the response variable is defined, while the right side
  outlines the fixed effect terms. Should only a single formula be
  provided - whether as a list with one item or as a standalone
  formula - it is inferred that a conventional univariate joint model is
  being constructed. Terms whose basis/contrasts depend on the data they
  are computed from – [`poly()`](https://rdrr.io/r/stats/poly.html) in
  its default orthogonal mode,
  [`splines::ns()`](https://rdrr.io/r/splines/ns.html)/
  [`splines::bs()`](https://rdrr.io/r/splines/bs.html), and
  [`factor()`](https://rdrr.io/r/base/factor.html) – are fine in a
  continuous biomarker's `long_sub_fixed` formula: the basis and factor
  levels the model was fit with are reused when estimating `Sigma_fit`
  and at prediction time. In `long_sub_random`, or in an ordinal
  biomarker's `long_sub_fixed`, they trigger a warning, because there
  the design matrix is still rebuilt from the data at hand (a single
  patient's rows at prediction time), so the basis can silently disagree
  with the one used to fit the model. Prefer `poly(..., raw = TRUE)`,
  `I(x^2)`, [`log()`](https://rdrr.io/r/base/Log.html),
  [`sqrt()`](https://rdrr.io/r/base/MathFun.html), or other terms that
  do not depend on the surrounding data there.

- long_sub_random:

  A list of one-sided formulas that define the model for the random
  effects of each longitudinal outcome. The number of items in this
  `list` should match the length of `formLongFixed`.

- biomarker_type:

  Optional character vector, one entry per longitudinal outcome, each
  either `"continuous"` or `"ordinal"`. When supplied, it always takes
  priority. When `NULL` (the default), each biomarker's type is
  auto-detected from its own response column in `data_fit_all`: a
  `factor`/`ordered factor` response is treated as `"ordinal"`, anything
  else as `"continuous"`. If every biomarker is continuous (the original
  use case), fitting proceeds exactly as before with no behavior change
  whatsoever. If at least one biomarker is ordinal, that biomarker is
  instead fit with
  [`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html) (a
  probit cumulative link mixed model, requiring the optional ordinal
  package) and folded into the shared random-effects covariance matrix
  via a Gaussian-copula extension of the EM algorithm – see
  [`longitudinalSubCopula()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSubCopula.md)
  for implementation details.

## Value

An object of class `"longitudinalSub.BJM"`, a named list with elements:

- lfit:

  A list of fitted univariate linear mixed models, one per longitudinal
  outcome, each obtained via
  [`lme`](https://rdrr.io/pkg/nlme/man/lme.html).

- Sigma_fit:

  The estimated variance-covariance matrix of the random effects in the
  multivariate linear mixed model.

- long_sub_fixed:

  The `long_sub_fixed` argument, as supplied.

- long_sub_random:

  The `long_sub_random` argument, as supplied.

- xlevels:

  A list, one element per longitudinal outcome, of the factor levels
  observed in the full training data for that outcome. Used internally
  at prediction time so that
  [`poly()`](https://rdrr.io/r/stats/poly.html)/
  [`splines::ns()`](https://rdrr.io/r/splines/ns.html)/[`splines::bs()`](https://rdrr.io/r/splines/bs.html)/[`factor()`](https://rdrr.io/r/base/factor.html)
  terms in `long_sub_fixed` reuse the basis/contrasts fit at training
  time instead of recomputing one from a small, patient-specific slice
  of data.

## Examples

``` r
# \donttest{

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

# Complete case analysis
data_fit_all = list()
for(i in seq_len(length(long_sub_fixed))){
  data_fit_all[[i]] = pbc3[pbc3$status3 == 1, ]
}

# fitting longitudinal submodel
long_fit_all = longitudinalSub(data_fit_all, long_sub_fixed, long_sub_random)

# poly() in its default orthogonal mode, splines::ns()/bs(), and
# factor() are safe in a continuous biomarker's long_sub_fixed: the
# terms/xlevels/contrasts fit on the full training data are cached and
# reused when estimating Sigma_fit and at prediction time, instead of
# being recomputed from each patient's small per-prediction slice.
long_fit_poly = longitudinalSub(
  pbc3[pbc3$status3 == 1, ],
  serBilir ~ year + poly(age, 2) + factor(sex) + years,
  ~ year | id)

# A few more nonlinear-term styles. None of these would trigger the
# warning even in long_sub_random, because their basis doesn't depend on the surrounding
# data at all -- there's simply nothing that could disagree between the
# full training data and the small per-patient slice used at prediction
# time.
data_fit_extra_terms = pbc3[pbc3$status3 == 1, ]
# A transformed-time column can also be built once with survivalTrans()
# and merged onto the fitting data ahead of time -- by the time
# long_sub_fixed sees it, it is just an ordinary numeric column, same as
# any of the other terms below.
data_fit_extra_terms$Tyears1 = survivalTrans(c(1, 3, 5, 7))$survival_trans_function[[1]](
  data_fit_extra_terms$year)

long_fit_nonlinear = list(
  longitudinalSub(data_fit_extra_terms, serBilir ~ year + I(year^2) + age + sex,
                   ~ year | id),
  longitudinalSub(data_fit_extra_terms, serBilir ~ year + I(year^2) * age + sex,
                   ~ year | id),
  longitudinalSub(data_fit_extra_terms, serBilir ~ year + log(year + 1) + age + sex,
                   ~ year | id),
  longitudinalSub(data_fit_extra_terms, serBilir ~ year + sqrt(year) + age + sex,
                   ~ year | id),
  longitudinalSub(data_fit_extra_terms, serBilir ~ poly(year, 2, raw = TRUE) + age + sex,
                   ~ year | id),
  longitudinalSub(data_fit_extra_terms, serBilir ~ year + Tyears1 + age + sex,
                   ~ year | id)
)

# }
```
