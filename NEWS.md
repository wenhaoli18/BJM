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
