# Convert an nlme-style random-effects formula to an lme4-style bar term

[`longitudinalSub()`](https://liwh0904.github.io/BJM/reference/longitudinalSub.md)'s
`long_sub_random` argument uses
[`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html)'s
`~ terms | group` formula convention.
[`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html) (used to
fit the marginal model for a categorical/ ordinal biomarker – see
[`fit_marginal_ordinal()`](https://liwh0904.github.io/BJM/reference/fit_marginal_ordinal.md))
instead expects random effects written as an `lme4`-style
`(terms | group)` bar term embedded directly in the model formula. This
helper translates one into the other so that the exact same
`long_sub_random` argument drives both
[`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) (continuous
biomarkers) and
[`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html) (ordinal
biomarkers), with no second random-effects argument for the user to keep
in sync.

## Usage

``` r
nlme_random_to_lme4_bars(random_formula)
```

## Value

A length-1 character string, e.g. `"(year | id)"`.
