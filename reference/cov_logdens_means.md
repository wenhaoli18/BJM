# Gaussian log-densities of many points under many means

For each mean vector in `means`, the log-density of every row of `x`,
with covariance factored by
[`cov_factor()`](https://wenhaoli18.github.io/BJM/reference/cov_factor.md).
Used by
[`conditionalYTBio()`](https://wenhaoli18.github.io/BJM/reference/conditionalYTBio.md)/[`conditionalYDTBio()`](https://wenhaoli18.github.io/BJM/reference/conditionalYDTBio.md),
where the rows of `x` are the candidate biomarker values and the means
are the survival-time grid points; the covariance does not depend on
either, so it is factored once per patient.

## Usage

``` r
cov_logdens_means(fac, x, means)
```

## Arguments

- fac:

  The result of
  [`cov_factor()`](https://wenhaoli18.github.io/BJM/reference/cov_factor.md).

- x:

  Matrix with one row per point, `N` columns.

- means:

  List of length-`N` mean vectors.

## Value

A list the same length as `means`; element `k` holds the log-densities
of the rows of `x` under `means[[k]]`.
