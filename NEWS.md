# BJM 0.3.0

## Breaking changes

* `dynamicPrediction()` has been renamed to `predictRisk()`, to read as a
  verb-first pair with the new `predictLongitudinal()` (below) rather than
  as one specifically-named function alongside a generically-named one.
  Existing code that calls `dynamicPrediction()` positionally or by name
  needs to switch to `predictRisk()`; the argument list, return value, and
  `"predictRisk.BJM"` (formerly `"dynamicPrediction.BJM"`) class are
  otherwise unchanged. `checkBandcountConvergence()`'s `predict_fun`
  argument accepts `predictRisk` in place of the old `dynamicPrediction`.

## New features

* New `imputeLongitudinal()` fills in missing longitudinal biomarker values
  before `longitudinalSub()` runs, instead of relying on complete-case
  analysis. `longitudinalSub()` drops rows missing a biomarker's own
  covariates, and then keeps only subjects who have at least one
  observation of *every* biomarker -- so with irregular/interrupted
  follow-up (missing at random given the observed covariates and
  biomarkers), a subject missing just one of several biomarkers at a visit,
  or missing a biomarker's measurements entirely, is dropped from every
  biomarker's fit, not just the one it is missing. `imputeLongitudinal()`
  fits a MIWAE (Missing data Importance-Weighted AutoEncoder; Mattei &
  Frellsen, 2019) jointly across the supplied biomarkers and draws one or
  more plausible completions via self-normalized importance resampling; the
  completed data can be passed straight into `longitudinalSub()` in place
  of the original `data_fit_all`. A second backend, selected with
  `method = "diffusion"`, fits a conditional denoising diffusion
  probabilistic model (DDPM; Ho, Jain & Abbeel, 2020) instead of a VAE,
  trained with the same observed-entries-only masking principle as the
  MIWAE objective (matching MissDiff, Ouyang et al., 2023), and imputes
  missing cells with a RePaint-style (Lugmayr et al., 2022) reverse-
  diffusion sampler that reproduces observed values exactly and only
  extrapolates the missing ones. Both backends share the same call shape,
  return shape, and `diagnostics` output, so switching between them is a
  one-argument change. Requires the optional `torch` package (listed in
  `Suggests`, not `Imports`, so installing/using BJM without it is
  unaffected); every call into `torch` is fully namespaced. The returned
  `diagnostics` also report a classical single-regression-imputation
  baseline alongside the deep-generative completions, so all of these can
  be compared on a given dataset -- deep generative imputation needs
  enough data to fit reliably, and is not automatically the better choice
  at every sample size.
* New `poolLongitudinalSub()` combines `longitudinalSub()` fits across the
  multiple completed datasets returned by `imputeLongitudinal(..., impute =
  "multiple")`, using Rubin's rules (Rubin, 1987) with the Barnard & Rubin
  (1999) small-sample degrees-of-freedom adjustment. `impute = "single"`
  (the default) fills `data_fit_all` with the across-draw mean of the
  generative model's completions and returns one completed dataset -- fast,
  but fitting `longitudinalSub()` on it treats every imputed value as if it
  had been observed, so the resulting standard errors do not reflect the
  uncertainty from not actually knowing the missing values. Fitting
  `longitudinalSub()` once per completion in `data_fit_all_list` (`impute =
  "multiple"`) and combining the fixed-effect estimates with
  `poolLongitudinalSub()` instead propagates that uncertainty into the
  pooled standard errors, degrees of freedom, and p-values, and reports a
  fraction-of-missing-information (`fmi`) per coefficient so it is visible
  which estimates were affected most by the missingness. `print()`-ing the
  result shows a coefficient table per biomarker alongside the FMI values.
* `longitudinalSub()` now accepts categorical (binary/ordinal) biomarkers
  alongside continuous ones, via a Gaussian-copula extension. Every
  biomarker's type -- `"continuous"` or `"ordinal"` -- is auto-detected from
  whether its response column in `data_fit_all` is a factor, or can be set
  explicitly with the new `biomarker_type` argument (which always takes
  priority over auto-detection when supplied). If every biomarker is
  continuous, fitting is completely unchanged: this is the same
  `nlme::lme()`-based code path as before, byte-for-byte. If at least one
  biomarker is ordinal, it is instead fit as a probit cumulative link mixed
  model via `ordinal::clmm()` (requiring the optional `ordinal` package,
  listed in `Suggests`, not `Imports`), representing each observed category
  as an interval on an underlying continuous latent Gaussian score --
  mathematically the same object as a Gaussian copula linking a discrete
  margin to the model's other (continuous or latent) margins. The shared
  random-effects covariance matrix across every biomarker, continuous and
  ordinal alike, is then re-estimated jointly by an ECM
  (Expectation-Conditional-Maximization) algorithm: one new "inner" E-step
  imputes each ordinal observation's latent score (its truncated-normal
  conditional mean given the fitted thresholds and current random-effect
  prediction) before the existing "outer" E-step/M-step -- unchanged from
  `longitudinalSub()`'s original EM -- re-estimates the covariance matrix.
  This inner-imputation step is a deterministic moment-matching plug-in
  (not a full Bayesian/MCEM draw), a documented v1 simplification.
  `predictRisk()` now also works directly on a fitted mixed-type
  model: its internal conditional-density calculation dispatches to a
  Gaussian-copula-aware variant (`conditionalYTCopula()`/
  `conditionalYDTCopula()`, internal, not exported) whenever
  `long_fit_all$biomarker_type` records at least one ordinal marker. Given
  a patient's observed continuous values and observed ordinal categories,
  this variant factors the joint density as the exact multivariate-normal
  density of the continuous sub-vector times the multivariate-normal box
  probability (`mvtnorm::pmvnorm()`) of the ordinal latent scores falling
  in their category-implied intervals, evaluated at the ordinal block's
  distribution conditional on the continuous observations (standard
  multivariate-normal conditioning). For an all-continuous fit
  (`biomarker_type` is `NULL`), dispatch falls through unchanged to the
  original `conditionalYT()`/`conditionalYDT()`, so no existing
  `predictRisk()` behavior is affected. `predictLongitudinal()`
  (predicting a future biomarker *value*, rather than event risk) now also
  dispatches to copula-aware variants the same way -- both its own
  conditional-density denominator (`conditionalYTCopula()`/
  `conditionalYDTCopula()`, shared with `predictRisk()`) and new
  `conditionalYTBioCopula()`/`conditionalYDTBioCopula()` numerators
  (internal, not exported) -- when predicting either a continuous **or an
  ordinal** biomarker from a mixed-type fit; the other jointly-fit
  biomarkers may freely be continuous or ordinal either way. Predicting the
  future *category* of an ordinal biomarker needs no separate code path:
  `conditionalYTBioCopula()`/`conditionalYDTBioCopula()` already bracket the
  candidate category's row between its cumulative-link thresholds via the
  same `biomarker_type`-driven machinery used for every other ordinal row in
  the joint density, and `predictLongitudinal()` represents an ordinal
  `bio_i`'s `Y_all`/`Y_predict` as integer category codes (`1:K`, in
  threshold order) rather than a numeric grid, with `Y_all` carrying a
  `"category_labels"` attribute giving the matching level names; `Y_density`
  is then that category's predicted probability rather than a density.
  `bandcount3` (the candidate-grid resolution) does not apply to an ordinal
  `bio_i`, whose grid is fixed at its category count, and is ignored in that
  case.
* New `predictLongitudinal()` is now the single entry point for predicting
  future biomarker value(s), replacing the two previously separate
  functions `dynamicPredictionBio()`/`dynamicPredictionBioAll()` (now
  internal, un-exported helpers behind it -- direct calls to either from
  existing scripts will need to switch to `predictLongitudinal()`). Pass a
  `bio_i` naming exactly **one** biomarker to predict just that one, which
  returns a single `dynamicPredictionBio.BJM` object; pass `bio_i` naming
  **more than one** biomarker, or leave it at the default `NULL` (meaning
  every biomarker in `long_fit_all`), to predict several at once, which
  instead returns a named list of such objects (one per biomarker, named by
  that biomarker's response variable, classed `dynamicPredictionBioAll.BJM`)
  and computes the bio_i-*independent* pipeline stages (restricting to
  at-risk patients, the survival-side integration grid, and the denominator
  conditional density) only **once**, reusing it across every requested
  biomarker instead of recomputing it once per biomarker -- redundant work
  that is especially costly under the Gaussian-copula path above, where
  that denominator involves `mvtnorm::pmvnorm()` Monte-Carlo evaluations.
  `bio_i` may name continuous and/or ordinal biomarkers, in any
  combination. `bandcount2` (the shared survival grid) is auto-tuned once,
  from a single representative biomarker, when left at its default
  `"auto"`, for the multi-biomarker case; `bandcount3` (each biomarker's
  own candidate-value grid) is then auto-tuned separately per biomarker,
  reusing the already-computed shared stage rather than rebuilding it on
  every doubling round -- except for any ordinal biomarker, whose candidate
  grid is fixed at its category count rather than controlled by
  `bandcount3`, so tuning is skipped for it and the returned `"bandcount3"`
  attribute records `NA` for that biomarker. `checkBandcountConvergence()`'s
  `predict_fun` argument now accepts `predictLongitudinal` (in place of the
  now-internal `dynamicPredictionBio`) for the single-biomarker case.

## Bug fixes

* `conditionalYT()` and `conditionalYDT()` -- the internal density functions
  behind `predictRisk()`'s risk-probability output whenever 2+
  longitudinal biomarkers are jointly fit -- computed the joint Gaussian
  density of the stacked observation vector across all markers via a "trace
  trick": reducing the quadratic form `(Y - mu)' Sigma_all_solve (Y - mu)`
  to `sum(diag(...))` of a smaller (number-of-markers-sized) matrix. That
  trace only sums the *diagonal* (same-marker) blocks of the reduced
  matrix, which silently discarded every *cross*-marker contribution of
  `Sigma_all_solve` -- i.e. it implicitly treated the biomarkers as
  conditionally independent given the random effects, even though
  `Sigma_all` is deliberately built from the fitted joint random-effects
  covariance (`Sigma_fit`) specifically to capture correlation *between*
  markers. Since `Sigma_fit` is generally not block-diagonal across
  markers, this under- or over-stated `predictRisk()`'s risk
  probabilities whenever the jointly-fit biomarkers had correlated random
  effects -- the ordinary case for a joint model, not an edge case.
  `dynamicPredictionBio()` was not affected: `conditionalYTBio()`/
  `conditionalYDTBio()` already computed this density directly via
  `mvtnorm::dmvnorm()` on the full covariance. Fixed by computing the
  quadratic form directly on the full stacked mean/observation vectors
  instead of via the trace reduction; verified against an independent
  `mvtnorm::dmvnorm()` reference on a real two-biomarker `pbc3` fit with
  correlated random effects (`test-conditional-density-correctness.R`).
  This changes the numeric value (but not the sign, scale of magnitude, or
  validity) of `predictRisk()`'s output for any existing multi-marker
  fit with correlated random effects; the golden-master regression
  baselines (`testdata/baseline.rds`, `testdata/baseline_noCR.rds`) have
  been regenerated to reflect the corrected values.

* `survivalSub()`'s internal `coxph()` call did not pass `x = TRUE, y =
  TRUE`, so the fitted model did not carry its own design matrix and
  response. Downstream functions that call `survival::basehaz()` or
  `survival::survfit()` on this fit (e.g. `predictRisk()` at
  prediction time) can, in that case, need to reconstruct the model frame
  by re-evaluating the fit's captured call -- but they do so in the
  environment of the fit's *formula*, i.e. the caller's environment, not
  `survivalSub()`'s own execution environment. Since `survivalSub()`
  passes its `data_survival_fitting` argument to `coxph()` by that same
  name, this happened to work whenever the caller's data object was also
  literally named `data_survival_fitting` (as in every example, vignette,
  and test shipped with the package) and failed with `object
  'data_survival_fitting' not found` for any other variable name. Fixed
  by passing `x = TRUE, y = TRUE` to the internal `coxph()` call, so the
  fit is self-contained and this reconstruction is never needed.

* `predictRisk()`/`dynamicPredictionBio()` failed with `"non-conformable
  arrays"` when predicting for a subject who had exactly one longitudinal
  observation to condition on, in a single-biomarker (univariate) model.
  The internal helpers `build_conditional_design()` (`R/conditionalDesign.R`)
  and `process_variance()` (`R/processVariance.R`) both built the residual
  covariance piece of the conditional variance as `diag(Sigma_vector)`,
  where `Sigma_vector` holds one residual variance per observed row for
  that subject. Base R's `diag()` is ambiguous on a length-1 numeric
  vector: instead of returning the intended 1x1 diagonal matrix, it
  interprets the single number as a *dimension* and returns an NxN
  identity matrix, which silently corrupted the dimensions of the
  downstream covariance matrix for exactly this one-observation case (any
  subject with 2+ observations was unaffected). Fixed by calling
  `diag(Sigma_vector, length(Sigma_vector))` in both places, which is
  unambiguous for every length, including 1.

# BJM 0.2.0

## Breaking changes

* All function arguments have been renamed to a single consistent
  `snake_case` convention. The package previously mixed `dot.case`,
  `camelCase`, and `PascalCase` for the same underlying concepts (and
  in a few cases used different names for the same concept in
  different functions). Code that calls these functions using
  positional arguments is unaffected; code that uses named arguments
  must update the names. The renames are:

  | Old name | New name |
  | --- | --- |
  | `data.survival.fitting` | `data_survival_fitting` |
  | `data.fit.all` | `data_fit_all` |
  | `data.plot.all` | `data_plot_all` |
  | `data.predict.all` | `data_predict_all` |
  | `data.predict.all.one` (`predictPlot`) | `data_predict_all_one` |
  | `data.predict.all.pre` (`riskPlot`) | `data_predict_all_pre` |
  | `prediction.time` | `prediction_time` |
  | `formMarginalSurv` | `form_marginal_surv` |
  | `formConditionalCR` | `form_conditional_cr` |
  | `survivalVariableAll` | `survival_variable_all` |
  | `survivalTransFunction` | `survival_trans_function` |
  | `LongSubFixed` | `long_sub_fixed` |
  | `LongSubRandom` | `long_sub_random` |

* Removed the `pbc2` example dataset. It was never referenced by any
  exported function, vignette, or test, and `pbc3` -- the dataset
  actually used throughout the package's examples and tests -- is not
  a duplicate of it (`pbc3` recodes `sex`, log-transforms several
  biomarkers, and adds the competing-risk/transformed-time columns
  `status3`, `status4`, `status5`, and `Tyears1`-`Tyears4`). Code that
  called `data(pbc2)` should switch to `data(pbc3)` and account for
  these differences.
* Removed the exported `print_survivalSub()`, `print_longitudinalSub()`,
  `print_dynamicPrediction()`, and `print_dynamicPredictionBio()`
  functions. Each only forwarded to the corresponding S3 `print.*.BJM`
  method, so calling `print()` on these objects (or letting them print
  at the console) is unaffected; only direct calls to the removed
  `print_*` functions need to switch to `print()`. `print_BJM()` has
  been renamed to `printBJM()`, since it prints a pair of independently
  fit sub-model objects together and does not participate in S3
  dispatch on a single `BJM`-classed object.

## New features

* `survivalSub()`, `longitudinalSub()`, `dynamicPrediction()`, and
  `dynamicPredictionBio()` now return named lists (e.g. `risk_prob_1`/
  `risk_prob_2`, `Y_predict`/`Y_density`/`Y_all`) with the fields
  documented in each function's `@return` block, instead of anonymous
  positional lists. Existing code indexing results with `x[[1]]`,
  `x[[2]]`, etc. continues to work unchanged.
* `survivalSub()`, `longitudinalSub()`, `dynamicPrediction()`, and
  `dynamicPredictionBio()` now validate their arguments up front and
  fail with a specific, actionable message (naming the offending
  argument) instead of a cryptic error from deep inside model-fitting
  or indexing code. `predictPlot()`, `riskPlot()`, and `cmtPlot()` got
  the same treatment.
* New `survivalTrans()` helper builds the `survival_variable_all`/
  `survival_trans_function` pair directly from a vector of cut points
  (e.g. `survivalTrans(c(1, 3, 5, 7))`), instead of requiring two
  hand-written, easy-to-misalign parallel lists.
* Every `data_*_all` argument across the pipeline (`data_fit_all`,
  `data_predict_all`, and now also `data_predict_all_one` in
  `predictPlot()` and `data_predict_all_pre` in `riskPlot()`) accepts a
  single bare `data.frame`, reused for every biomarker, instead of a
  repeated list, when all biomarkers share the same measurement data.
* `longitudinalSub()` now warns when a `long_sub_fixed` formula contains
  `poly()` (in its default orthogonal mode), `splines::ns()`,
  `splines::bs()`, or `factor()`. These terms recompute their
  basis/contrasts from whatever data they are given, but
  `dynamicPrediction()`/`dynamicPredictionBio()` rebuild the design
  matrix from a small, patient-specific slice of data at every point on
  the internal prediction grid, not the data the model was fit on -- so
  the basis silently disagrees with the one used at fitting time
  (wrong predictions) or `model.matrix()` fails outright when there are
  too few distinct values. Use `poly(..., raw = TRUE)`, `I(x^2)`,
  `log()`, `sqrt()`, or other terms that do not depend on the
  surrounding data instead.
* `bandcount1`/`bandcount2`/`bandcount3` (in `dynamicPrediction()`,
  `dynamicPredictionBio()`, `predictPlot()`, and `riskPlot()`) no longer
  need to be chosen by hand: they now default to `"auto"` instead of a
  fixed number. Under `"auto"`, the value is started small and doubled,
  comparing the returned predictions to the previous round, until the
  largest relative change drops below 1%, or 2 doublings have been
  tried (so resolving a bandcount costs at most 3 prediction calls, not
  an open-ended loop). `predictPlot()`/`riskPlot()` resolve their
  `"auto"` bandcount(s) once, using a representative probe call, rather
  than repeating the search on every point in their internal
  `horizon`/landmark-time loop. If a bandcount has still not converged
  after hitting this cap, a warning reports it and the result at the
  largest value tried is returned anyway (not an error). Pass an
  explicit number, as in previous package versions, to skip auto-tuning
  and use a fixed value instead.
* New `checkBandcountConvergence()` helper gives direct, manual control
  over the same doubling check that now runs automatically by default
  (e.g. to use a tighter tolerance, or more doublings, than the
  built-in `"auto"` search): it runs `dynamicPrediction()`/
  `dynamicPredictionBio()` once at the bandcount value(s) you supply and
  once more with those value(s) scaled up (by default, doubled), and
  reports the largest relative change in the returned predictions -- at
  the cost of exactly one extra prediction call.

## Bug fixes

* `riskPlot()` built its internal `data_predict_all` accumulator
  without initializing it first, so it could silently pick up a
  leftover object of the same name from the caller's environment.
  Fixed to initialize it explicitly, matching `predictPlot()`.
* The "wrap a single `data.frame` as a list" convenience documented for
  `longitudinalSub()`, `conditionalYT()`, `conditionalYTBio()`,
  `conditionalYDT()`, `conditionalYDTBio()`, and `process_variance()`
  never actually triggered, because `is.list()` is `TRUE` for
  data frames in R. Fixed the guard in all six places to also check
  `is.data.frame()`.
* `cmtPlot()`'s `id_variable` argument was dead code: three internal
  deduplication steps always looked up a literal column named
  `"id_variable"` instead of the column named by the argument, so a
  custom `id_variable` silently had no effect. Fixed to look up the
  specified column.
* `cmtPlot(condi_time2event = NULL)` crashed with `object 'plot_data'
  not found`, because the fallback that picks the midpoint of
  `time_variable` referenced an undefined variable instead of the
  actual `data_plot_all` argument. Fixed, so `condi_time2event = NULL`
  works as documented.
* `dynamicPrediction()`, `dynamicPredictionBio()`, `predictPlot()`, and
  `riskPlot()` only checked that each element of
  `survival_trans_function` was a function, never that it actually
  worked. A transform that throws an error, or returns a character
  value, a vector of the wrong length, or a non-finite value (e.g.
  `log(x)` at `x <= 0`), previously only surfaced as a cryptic failure
  deep inside the per-patient prediction grid, or, in the non-finite
  case, as silently corrupted predictions with no error at all. Fixed
  by test-calling every transform once, up front, on the supplied
  `prediction_time` and validating its output, before any of the
  (potentially expensive) prediction machinery runs.

## Internal changes

* `longitudinalSub()`'s documentation now includes a worked
  `poly()`/`splines::ns()`/`factor()` example, matched by an equivalent
  example in the package's example scripts, showing that the warning
  described above is a prompt to double-check the fitted basis, not a
  sign that predictions from these terms are wrong.
* Fixed `pbc3`'s documentation, which had copied `@usage`/`@format`
  tags from `pbc2` (`data(pbc2)`, "20 variables") instead of describing
  `pbc3` itself (`data(pbc3)`, 27 variables); the variable-by-variable
  `\describe` list was unaffected and already documented all 27 `pbc3`
  columns correctly.
* Added a `tests/testthat` suite, including golden-master
  characterization tests captured from the package's own `pbc3`
  example data, covering both the competing-risk and
  no-competing-risk code paths.
* Extracted shared per-patient design-matrix construction out of
  `conditionalYT()`/`conditionalYDT()` and
  `conditionalYTBio()`/`conditionalYDTBio()` into internal helpers in
  `R/conditionalDesign.R`.
* Extracted shared at-risk subsetting, integration-grid setup, and
  risk-probability clamping logic out of `dynamicPrediction()` and
  `dynamicPredictionBio()` into internal helpers in
  `R/dynamicPredictionShared.R`.
* No numeric behavior change is intended by the internal changes in
  this release; all changes were verified against golden-master
  baselines captured before refactoring.

# BJM 0.1.0

* Initial CRAN release.
