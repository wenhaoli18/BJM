# Call a single survival_trans_function element and validate its output

[`assert_survival_trans()`](https://liwh0904.github.io/BJM/reference/assert_survival_trans.md)
only probes each transform once, at a single representative time value,
before the prediction grid runs. That catches a transform that is broken
everywhere (wrong return type/length, or throws), but not one that only
misbehaves away from the probe point – e.g. `log(x - 10)`, which is fine
near the probe value but returns `NaN` once the internal prediction grid
(which can range up to `2 * max(observed survival time)`) reaches
`x <= 10`. This helper wraps every actual call site inside
[`conditionalYT()`](https://liwh0904.github.io/BJM/reference/conditionalYT.md)/
[`conditionalYTBio()`](https://liwh0904.github.io/BJM/reference/conditionalYTBio.md)/[`conditionalYDT()`](https://liwh0904.github.io/BJM/reference/conditionalYDT.md)/[`conditionalYDTBio()`](https://liwh0904.github.io/BJM/reference/conditionalYDTBio.md)
so a bad value is caught immediately, with a clear error, instead of
silently corrupting a data.frame column or surfacing later as a cryptic
"replacement has ... rows, data has ..." error.

## Usage

``` r
apply_survival_trans(fun, x, surv_i)
```
