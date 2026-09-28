# Assert that survival_variable_all/survival_trans_function are consistent

Shared input-validation helper for
[`predictRisk()`](https://liwh0904.github.io/BJM/reference/predictRisk.md),
[`dynamicPredictionBio()`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBio.md),
[`predictPlot()`](https://liwh0904.github.io/BJM/reference/predictPlot.md),
and
[`riskPlot()`](https://liwh0904.github.io/BJM/reference/riskPlot.md):
the two arguments must have matching length, and every transform must be
a function. When `probe_value` is supplied, every transform is also
test-called once on it, and must return a single, finite, non-missing
numeric value. Without this, a transform that throws an error, or
returns a character value, a length != 1 vector, or a non-finite value
(e.g. `log(x)` evaluated at `x <= 0`), would only surface deep inside
the per-patient prediction grid built by
[`conditionalYT()`](https://liwh0904.github.io/BJM/reference/conditionalYT.md)/
[`conditionalYDT()`](https://liwh0904.github.io/BJM/reference/conditionalYDT.md)/[`conditionalYTBio()`](https://liwh0904.github.io/BJM/reference/conditionalYTBio.md)/[`conditionalYDTBio()`](https://liwh0904.github.io/BJM/reference/conditionalYDTBio.md)
– as a cryptic error, or, worse, as silently corrupted data with no
error at all. The probe is a single call per transform, so it is cheap
even though the same transform is later called many times inside the
prediction grid.

## Usage

``` r
assert_survival_trans(
  survival_variable_all,
  survival_trans_function,
  probe_value = NULL
)
```
