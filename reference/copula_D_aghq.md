# Copula EM with adaptive Gauss-Hermite quadrature

Given the ordinal random effects \\b_o\\, the continuous ones \\b_c\\
have a Gaussian posterior with mean \\a + B b_o\\ and covariance
\\A_c\\, and the continuous data contribute a Gaussian factor in
\\b_o\\; only \\b_o\\ is integrated numerically. For each subject the
posterior mode of \\b_o\\ (Newton, warm-started from the previous
iteration) and the curvature there place a product Gauss-Hermite grid,
which gives \\E\[b b^T \mid y\]\\ block by block; their average over
subjects is the next \\D\\. Uses 7 nodes per dimension up to two ordinal
random effects and 5 beyond. The EM is accelerated by SQUAREM (Varadhan
and Roland, 2008), falling back to plain EM steps where an extrapolated
\\D\\ is not positive definite; convergence is a change in \\D\\, scaled
by \\\sqrt{D\_{jj} D\_{kk}}\\ so that it does not depend on the
biomarkers' units, below `tol.em` over one accelerated cycle (three EM
updates).

## Usage

``` r
copula_D_aghq(D, setup, tol.em = 1e-04, max.iter = 500, verbose = FALSE)
```

## Arguments

- D:

  Starting covariance.

- setup:

  Output of
  [`copula_estep_setup()`](https://wenhaoli18.github.io/BJM/reference/copula_estep_setup.md).

- tol.em, max.iter:

  Convergence tolerance on `D` and iteration cap of the quadrature EM.

- verbose:

  Print the change in `D` at every iteration.

## Value

A list with `D`.
