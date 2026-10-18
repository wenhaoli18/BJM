# `mixed_density_prob_copula()` for many means under one covariance

Factors the parts of
[`mixed_density_prob_copula()`](https://wenhaoli18.github.io/BJM/reference/mixed_density_prob_copula.md)
that depend only on the covariance – the inverse and log-determinant of
the continuous block, and the conditional covariance of the ordinal
block – once, and returns a function that evaluates the log density for
many means, observed values and threshold bounds at once (one per
column), e.g. a patient's density at every integration grid point.

The ordinal box probability is computed by adaptive Gauss-Hermite
quadrature
([`box_prob_quadrature()`](https://wenhaoli18.github.io/BJM/reference/box_prob_quadrature.md))
when the patient's ordinal rows involve a single ordinal random effect
(e.g. one ordinal biomarker with a random intercept). Given the
continuous values, the ordinal latent scores then have covariance \\l
l^T + I\\, \\l l^T\\ from the random effect and \\I\\ from the latent
residuals, so they are independent given the random effect, and the
probability is a one-dimensional integral of a product of univariate
normal probabilities. It used to be computed by
[`mvtnorm::pmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/pmvnorm.html)'s
Monte Carlo algorithm, which dominated the run time of every prediction
with an ordinal biomarker, differed between runs (relative error around
\\10^{-4}\\), and needed
[`set.seed()`](https://rdrr.io/r/base/Random.html) to be reproducible.
With two or more ordinal random effects `pmvnorm()` is still used: a
product quadrature grid then needs so many nodes to match its accuracy
(a random intercept with a large variance makes the integrand lopsided)
that it is no faster. It is also used when the covariance does not have
this structure, and when `box = "pmvnorm"`.

## Usage

``` r
mixed_density_copula_prep(
  Sigma_all,
  cc_idx,
  oo_idx,
  box = c("quadrature", "pmvnorm")
)
```

## Arguments

- Sigma_all:

  Full stacked covariance matrix for this patient (as returned by
  [`build_conditional_design_copula()`](https://wenhaoli18.github.io/BJM/reference/build_conditional_design_copula.md)).

- cc_idx, oo_idx:

  Integer index vectors (into the stacked row order) of the continuous-
  and ordinal-marker rows respectively.

- box:

  `"quadrature"` (default) or `"pmvnorm"`.

## Value

A function `(mu_full, y_all_vec, alpha_lower_vec, alpha_upper_vec)`
returning the log densities. `mu_full` is a vector or a matrix with one
column per evaluation; the other three may be vectors (the same for
every column) or matrices of the same shape.
