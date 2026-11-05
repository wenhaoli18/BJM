# CRAN submission comments

This is a resubmission fixing the issue reported on the CRAN check page
for BJM 0.3.0 (to be fixed before 2026-11-23): under the MKL additional
check, two `expect_identical()` comparisons in
`tests/testthat/test-longitudinalSubCopula.R` failed because refitting
the same `nlme::lme()` model gave results differing in about the 12th
significant digit. The tests that compare two fits of the same model now
use `expect_equal()` (testthat's default numerical tolerance) instead.
There is no change to the package code. I could not run an MKL build
locally.

## R CMD check results

win-builder, R-devel (2026-11-03 r90835 ucrt), Windows:

0 errors | 0 warnings | 1 note

* `checking CRAN incoming feasibility ... NOTE` -- "Days since last
  update: 4". This update fixes the issue reported on the CRAN check page
  shortly after 0.3.0 was published.

`R CMD check --as-cran` (R 4.4.1, macOS 14.4 aarch64, local) reports:

0 errors | 0 warnings | 3 notes

* `checking CRAN incoming feasibility ... NOTE` -- "Days since last
  update", as above.

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

`testthat` (3rd edition), as run by `R CMD check`: 801 passed, 0 failed,
0 warnings, 59 skipped on CRAN, in about 1 minute locally. The skipped
tests are long-running ones (reference results, plots, and repeated
fits of ordinal biomarkers with `ordinal::clmm()`), skipped to keep the
check time down; each model path still runs at least one test on CRAN,
and all 1098 tests pass with `NOT_CRAN=true`.

## Downstream dependencies

This is a leaf package with no known reverse dependencies on CRAN.
