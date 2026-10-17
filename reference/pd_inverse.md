# Inverse of a covariance matrix that may be close to singular

By the eigendecomposition, with eigenvalues floored at `1e-10` times the
largest, so a covariance estimate approaching the boundary (e.g. a
random-effect variance tending to 0) still gives a finite,
positive-definite precision.

## Usage

``` r
pd_inverse(V)
```

## Arguments

- V:

  Symmetric matrix.

## Value

Its (regularized) inverse.
