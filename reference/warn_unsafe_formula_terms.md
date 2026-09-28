# Warn about formula terms whose basis is recomputed from whatever data they are given

[`poly()`](https://rdrr.io/r/stats/poly.html) (in its default orthogonal
mode),
[`splines::ns()`](https://rdrr.io/r/splines/ns.html)/[`splines::bs()`](https://rdrr.io/r/splines/bs.html),
and [`factor()`](https://rdrr.io/r/base/factor.html) compute their
basis/contrasts from whatever data is passed to
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html). BJM's
prediction functions rebuild the design matrix from a small,
patient-specific slice of data at every point on the internal prediction
grid, which is not the data the model was fit on, so the basis
recomputed at prediction time silently does not match the one used at
fitting time (or, with too few distinct values,
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) fails
outright). `poly(..., raw = TRUE)`, `I(x^2)`,
[`log()`](https://rdrr.io/r/base/Log.html),
[`sqrt()`](https://rdrr.io/r/base/MathFun.html), and similar terms that
do not depend on the surrounding data are unaffected and are not
flagged.

## Usage

``` r
warn_unsafe_formula_terms(formula_list, arg_name)
```
