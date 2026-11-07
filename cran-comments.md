# CRAN submission comments

## Resubmission

Thank you for checking the previous submission of 0.3.1 under MKL. It
still failed there: switching from `expect_identical()` to
`expect_equal()` was not enough, because the fits being compared contain
`nlme::lme()`'s `apVar`, a finite-difference approximation to the
covariance of the variance parameters. Last-digit differences in the fit
change it by tens of percent. I also found that the second biomarker's
`optim()` fit can stop at relative differences of a few 1e-8, just above
`expect_equal()`'s default tolerance.

In this version, the three comparisons in
`tests/testthat/test-longitudinalSubCopula.R` drop `apVar` (which BJM does
not use) and use a tolerance of 1e-6. I could not run an MKL build
locally, so I reproduced the problem instead by refitting on data
perturbed in the 10th to 14th significant digit: `apVar` changed by up to
57% and the previous test failed, while the new comparison passed in all
90 such refits. It still fails for a genuinely different fit (other fixed
effects, other random effects, or data changed by 1e-4). There is no
change to the package code.

## Original submission

This update fixes the issue reported on the CRAN check page for BJM 0.3.0
(to be fixed before 2026-11-23): under the MKL additional check, tests
that refit the same `nlme::lme()` model twice and required the two fits to
be identical failed.

## R CMD check results

TODO before submitting: rerun win-builder on the final tarball and fill
in the result below.

win-builder, R-devel (TODO: date and revision), Windows:

0 errors | 0 warnings | 1 note

* `checking CRAN incoming feasibility ... NOTE` -- "Days since last
  update". This update fixes the issue reported on the CRAN check page
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
