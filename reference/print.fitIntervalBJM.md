# Print a fitIntervalBJM fit

A short summary of a
[`fitIntervalBJM`](https://wenhaoli18.github.io/BJM/reference/fitIntervalBJM.md)
fit: the imputation settings, the survival sub-model's hazard ratios and
the longitudinal sub-model's fixed effects (pooled over the imputations,
with the fraction of missing information). The components themselves are
in the returned list (e.g. `x$pooled`, `x$trace`).

## Usage

``` r
# S3 method for class 'fitIntervalBJM'
print(x, digits = 4, ...)
```

## Arguments

- x:

  A `"fitIntervalBJM"` object.

- digits:

  Number of significant digits.

- ...:

  Not used.

## Value

`x`, invisibly.
