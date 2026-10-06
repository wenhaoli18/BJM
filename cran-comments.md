# CRAN submission comments

This is an update to the previously published version (0.2.0 on CRAN)
to version 0.3.0. See NEWS.md for the full list of changes, including
the rename of `dynamicPrediction()` to `predictRisk()`, the new
`predictLongitudinal()`, `imputeLongitudinal()` (optional `torch`
backend, in `Suggests`), and `poolLongitudinalSub()` functions,
support for ordinal biomarkers, and several bug fixes.

## Update shortly after the previous release

This submission follows 0.2.0 closely. It fixes a correctness issue in
the dynamic prediction functions: longitudinal measurements recorded
after `prediction_time` were not removed from the prediction data, so a
prediction could condition on future information. They are now dropped
(with a warning) before predicting.

## R CMD check results

`R CMD check --as-cran` (R 4.4.1, macOS, local) reports:

0 errors | 0 warnings | 2 notes

* `checking CRAN incoming feasibility ... NOTE` -- "Days since last
  update: 5"; see "Update shortly after the previous release" above.
* `checking for future file timestamps ... NOTE` -- "unable to verify
  current time"; a network/clock-verification check with no bearing on
  the package itself.

`testthat` (3rd edition), as run by `R CMD check`: 403 passed, 0 failed,
2 skipped on CRAN (long-running baseline-characterization tests).

## Downstream dependencies

This is a leaf package with no known reverse dependencies on CRAN.
