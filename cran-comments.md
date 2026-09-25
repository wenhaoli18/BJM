# CRAN submission comments

This is an update to the previously published version (0.1.0 on CRAN)
to version 0.2.0. See NEWS.md for the full list of changes, including
breaking argument-name and return-value changes, new features
(`survivalTrans()`, auto-tuned bandcounts, `checkBandcountConvergence()`,
bare-`data.frame` support), several bug fixes, and removal of the
unused `pbc2` example dataset (superseded by `pbc3`, which the
package's examples and tests already use throughout).

## R CMD check results

`R CMD check --as-cran` (R 4.4.1, macOS, local) reports:

0 errors | 0 warnings | 2 notes

* `checking for future file timestamps ... NOTE` -- "unable to verify
  current time"; a network/clock-verification check with no bearing on
  the package itself.
* `checking HTML version of manual ... NOTE` -- Tidy HTML5 validation
  warnings (e.g. `<main>` not recognized, "onload" flagged as a
  proprietary attribute, `<table>` lacking a `summary` attribute) in
  Rd-to-HTML rendering. These come from R's own Rd-to-HTML converter
  output, not from any HTML authored in this package, and do not
  affect the rendered manual's usability.

`devtools::test()` (testthat, 3rd edition): 128 passed, 0 failed.

## Downstream dependencies

This is a leaf package with no known reverse dependencies on CRAN.
