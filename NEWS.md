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

## Internal changes

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
