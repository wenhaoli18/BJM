# Gaussian log-density from a factored covariance

Evaluates the `N`-variate normal log-density with covariance factored by
[`cov_factor()`](https://wenhaoli18.github.io/BJM/reference/cov_factor.md)
at each column of a residual matrix (observation minus mean).

## Usage

``` r
cov_logdens(fac, E)
```

## Arguments

- fac:

  The result of
  [`cov_factor()`](https://wenhaoli18.github.io/BJM/reference/cov_factor.md).

- E:

  Residual matrix, `N x m` (or a length-`N` vector).

## Value

A length-`m` vector of log-densities. If the covariance is not positive
definite (`"dense"` only), `-Inf`, or `Inf` for a zero residual, as
[`mvtnorm::dmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/Mvnorm.html)
returns.
