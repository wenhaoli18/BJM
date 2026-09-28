# Assert that a formula's variables are all present in a data.frame

Shared input-validation helper: catches missing columns before they
surface as an opaque error from deep inside
[`coxph()`](https://rdrr.io/pkg/survival/man/coxph.html), `lme()`, or
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html).

## Usage

``` r
assert_vars_in_data(vars, data, source_name, data_name)
```
