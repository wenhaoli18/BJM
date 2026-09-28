# Assert that a bandcount argument is a positive number or "auto"

Shared input-validation helper for the `bandcount1`/
`bandcount2`/`bandcount3` arguments of
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md),
[`dynamicPredictionBio()`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md),
[`predictPlot()`](https://wenhaoli18.github.io/BJM/reference/predictPlot.md),
and
[`riskPlot()`](https://wenhaoli18.github.io/BJM/reference/riskPlot.md),
which now accept either an explicit positive number (the original
behavior) or the literal string `"auto"` to have the value chosen
automatically (see
[`auto_tune_bandcount()`](https://wenhaoli18.github.io/BJM/reference/auto_tune_bandcount.md)).

## Usage

``` r
assert_bandcount(x, arg_name)
```
