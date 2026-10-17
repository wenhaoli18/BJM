# Formula terms whose values depend on the other rows of the data

The calls in `expr` to `funs` – by default
[`poly()`](https://rdrr.io/r/stats/poly.html) (unless `raw = TRUE`),
[`splines::ns()`](https://rdrr.io/r/splines/ns.html),
[`splines::bs()`](https://rdrr.io/r/splines/bs.html) and
[`factor()`](https://rdrr.io/r/base/factor.html) – as deparsed strings.
See
[`warn_unsafe_formula_terms()`](https://wenhaoli18.github.io/BJM/reference/warn_unsafe_formula_terms.md).

## Usage

``` r
unsafe_formula_calls(expr, funs = c("poly", "ns", "bs", "factor"))
```

## Value

A character vector, without duplicates.
