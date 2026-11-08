# Warn once when the longitudinal sub-model ignores many censored subjects

Called by
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md):
if the longitudinal fit has no more subjects than the survival fit has
events while at least 30\\ the survival subjects are censored, it was
most likely fit on the subjects with an event only, and risks can be
overestimated (see
[`fitImputedBJM`](https://wenhaoli18.github.io/BJM/reference/fitImputedBJM.md)).
Warns once per session; set `options(BJM.check_complete_case = FALSE)`
to silence it.

## Usage

``` r
check_complete_case_fit(long_fit_all, survival_fit_all)
```
