# Joint random-effects covariance with ordinal biomarkers (Gaussian copula)

Re-estimates the shared random-effects covariance `D` across all
biomarkers by maximum likelihood, holding each biomarker's fixed
effects, residual variance and (for an ordinal biomarker) thresholds at
their separate-fit values, as
[`longitudinalSubVar()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSubVar.md)
does for continuous ones. An ordinal observation is the category of a
latent score \\Z = \eta + \epsilon\\, \\\epsilon \sim N(0, 1)\\ (probit
link, Albert & Chib), observed only through the thresholds that bracket
it, so the posterior of the random effects given the data is not
Gaussian and the EM's \\E\[b b^T \mid y\]\\ has no closed form. It is
computed by

- [`copula_D_aghq()`](https://wenhaoli18.github.io/BJM/reference/copula_D_aghq.md),
  an EM with adaptive Gauss-Hermite quadrature over the ordinal
  biomarkers' random effects (the continuous biomarkers' random effects
  are integrated out exactly), when those number at most
  `max_quadrature_dim`; deterministic, with the usual convergence check;

- otherwise
  [`copula_D_saem()`](https://wenhaoli18.github.io/BJM/reference/copula_D_saem.md),
  a stochastic-approximation EM (SAEM) that draws the latent scores,
  whose cost does not grow with that dimension.

Both target the same maximum-likelihood estimate. The earlier version
plugged in the mean of each latent score given its category, which
ignored the score's own variance and underestimated the ordinal
biomarkers' random-effects variances (by about 40\\ for a random slope
in simulations).

## Usage

``` r
longitudinalSubVarCopula(
  thetaLong,
  l,
  biomarker_type,
  thresholds,
  tol.em = 1e-04,
  max.iter = 500,
  max_quadrature_dim = 3,
  verbose = FALSE
)
```

## Arguments

- thetaLong:

  List with the stacked fixed effects `beta`, the starting `D` and the
  residual variances `sigma2` (1 for an ordinal biomarker).

- l:

  Per-subject design bookkeeping built by
  [`longitudinalSubCopula()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSubCopula.md);
  `l$yik[[m]]` holds an ordinal biomarker's category codes.

- biomarker_type, thresholds:

  As stored in the fit.

- tol.em, max.iter:

  Convergence tolerance on `D` and iteration cap of the quadrature EM.

- max_quadrature_dim:

  Largest number of ordinal random effects handled by quadrature.

- verbose:

  Print the change in `D` at every iteration.

## Value

A list with the estimate `D`.
