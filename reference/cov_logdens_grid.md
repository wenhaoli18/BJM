# Log densities of every candidate outcome vector at every grid point

For
[`conditionalYTBio()`](https://wenhaoli18.github.io/BJM/reference/conditionalYTBio.md)/[`conditionalYDTBio()`](https://wenhaoli18.github.io/BJM/reference/conditionalYDTBio.md):
the log density of each row of `x` (one per candidate value of the
predicted biomarker) with mean `mu[, k]`, for every grid point `k`, in
one
[`cov_logdens()`](https://wenhaoli18.github.io/BJM/reference/cov_logdens.md)
call.

## Usage

``` r
cov_logdens_grid(fac, x, mu)
```

## Arguments

- fac:

  The patient's factored covariance
  ([`cov_factor()`](https://wenhaoli18.github.io/BJM/reference/cov_factor.md)).

- x:

  Matrix, one row per candidate value, one column per stacked
  observation.

- mu:

  Matrix of means, one row per stacked observation, one column per grid
  point.

## Value

A `nrow(x)` x `ncol(mu)` matrix.
