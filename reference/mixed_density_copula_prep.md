# `mixed_density_prob_copula()` for many means under one covariance

Factors the parts of
[`mixed_density_prob_copula()`](https://wenhaoli18.github.io/BJM/reference/mixed_density_prob_copula.md)
that depend only on the covariance – the inverse and log-determinant of
the continuous block, and the conditional covariance of the ordinal
block – once, and returns a function of the mean, observed values and
threshold bounds that does the rest with the same arithmetic. Used to
evaluate a patient's density at every integration grid point.

## Usage

``` r
mixed_density_copula_prep(Sigma_all, cc_idx, oo_idx)
```

## Arguments

- Sigma_all:

  Full stacked covariance matrix for this patient (as returned by
  [`build_conditional_design_copula()`](https://wenhaoli18.github.io/BJM/reference/build_conditional_design_copula.md)).

- cc_idx, oo_idx:

  Integer index vectors (into the stacked row order) of the continuous-
  and ordinal-marker rows respectively.

## Value

A function `(mu_full, y_all_vec, alpha_lower_vec, alpha_upper_vec)`
returning the log density.
