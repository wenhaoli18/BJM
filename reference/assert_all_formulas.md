# Assert that every element of a list is a formula

Shared input-validation helper for `long_sub_fixed`/
`long_sub_random`-style arguments, after they have been normalized to a
list (a bare formula is wrapped in
[`list()`](https://rdrr.io/r/base/list.html) by the caller before
calling this).

## Usage

``` r
assert_all_formulas(x, arg_name)
```
