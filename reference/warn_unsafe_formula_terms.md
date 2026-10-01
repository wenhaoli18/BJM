# Warn about formula terms whose basis is recomputed from whatever data they are given

[`poly()`](https://rdrr.io/r/stats/poly.html) (in its default orthogonal
mode),
[`splines::ns()`](https://rdrr.io/r/splines/ns.html)/[`splines::bs()`](https://rdrr.io/r/splines/bs.html),
and [`factor()`](https://rdrr.io/r/base/factor.html) compute their
basis/contrasts from whatever data is passed to
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html). For a
continuous biomarker's `long_sub_fixed` formula this is handled: the
basis `lme()` was fit with (its `terms`' `predvars`) and the training
factor levels are reused both when estimating `Sigma_fit` and at
prediction time, so such terms are not flagged there. They are still
rebuilt from the data at hand – a single patient's rows at prediction
time, the retained subjects when estimating `Sigma_fit` – for
`long_sub_random` formulas and for an ordinal biomarker's
`long_sub_fixed` formula, which is where this warns.
`poly(..., raw = TRUE)`, `I(x^2)`,
[`log()`](https://rdrr.io/r/base/Log.html),
[`sqrt()`](https://rdrr.io/r/base/MathFun.html), and similar terms that
do not depend on the surrounding data are never flagged.

## Usage

``` r
warn_unsafe_formula_terms(formula_list, arg_name)
```
