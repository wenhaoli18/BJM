# Assert that an object is a single positive integer

Shared input-validation helper for `n_cores`: catches a non-integer,
zero, negative, or non-scalar value before it reaches
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)'s own
(less informative) `mc.cores` validation.

## Usage

``` r
assert_positive_integer(x, arg_name)
```
