# CRAN submission comments

This is an update to the previously published version (0.2.0 on CRAN,
published 2026-09-25) to version 0.3.0. See NEWS.md for the full list of
changes, including the rename of `dynamicPrediction()` to `predictRisk()`,
the new `predictLongitudinal()`, `imputeLongitudinal()` (optional `torch`
backend, in `Suggests`), `poolLongitudinalSub()` and `fitIntervalBJM()`
functions, support for ordinal biomarkers and interval-censored event
times, and several bug fixes -- among them a correctness fix in the
dynamic prediction functions: longitudinal measurements recorded after
`prediction_time` were not removed from the prediction data, so a
prediction could condition on future information. They are now dropped
(with a warning) before predicting.

## R CMD check results

win-builder, R-devel (2026-10-14 r90673 ucrt), Windows:

0 errors | 0 warnings | 1 note

* `checking CRAN incoming feasibility ... NOTE` -- "Possibly misspelled
  words in DESCRIPTION: Shen, al, et". These are false positives from the
  method references in the Description field, written in the requested
  "authors (year) <doi:...>" form: "Shen" is an author's surname (Shen and
  Li (2021)), and "et al." abbreviates the author list of Li et al. (2026).

`R CMD check --as-cran` (R 4.4.1, macOS 14.4 aarch64, local) reports:

0 errors | 0 warnings | 2 notes

* `checking for future file timestamps ... NOTE` -- "unable to verify
  current time"; a network/clock-verification check with no bearing on
  the package itself.
* `checking HTML version of manual ... NOTE` -- validation messages such
  as "`<main>` is not recognized!" from the local HTML Tidy, which is the
  2006 build shipped with macOS and predates HTML5; not a problem in the
  package's Rd files.

Examples, including those in `\donttest{}`, run without error. The `torch`
examples run only when both the `torch` package and its LibTorch backend
are installed (`torch::torch_is_installed()`).

`testthat` (3rd edition), as run by `R CMD check`: 586 passed, 0 failed,
0 warnings, 53 skipped on CRAN, in about 1 minute locally. The skipped
tests are long-running ones (reference results, plots, and repeated
fits of ordinal biomarkers with `ordinal::clmm()`), skipped to keep the
check time down; each model path still runs at least one test on CRAN,
and all 854 tests pass with `NOT_CRAN=true`.

## Downstream dependencies

This is a leaf package with no known reverse dependencies on CRAN.
