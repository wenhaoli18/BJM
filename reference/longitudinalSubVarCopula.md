# EM re-estimation of the shared random-effects covariance (Gaussian-copula extension)

Extends
[`longitudinalSubVar()`](https://liwh0904.github.io/BJM/reference/longitudinalSubVar.md)'s
EM algorithm with one additional "inner" E-step
([`impute_latent_ordinal()`](https://liwh0904.github.io/BJM/reference/impute_latent_ordinal.md)),
run before the existing "outer" E-step (posterior mean/second-moment of
the random effects) at every iteration – an ECM
(Expectation-Conditional-Maximization) scheme. The outer E-step and
M-step formulas are otherwise completely unchanged from
[`longitudinalSubVar()`](https://liwh0904.github.io/BJM/reference/longitudinalSubVar.md):
only the stacked outcome vector `yi` fed into them differs, since
ordinal biomarkers' observed category codes are replaced by the inner
E-step's imputed latent scores before every outer-step update. As in
[`longitudinalSubVar()`](https://liwh0904.github.io/BJM/reference/longitudinalSubVar.md),
`sigma2` is never updated (for ordinal biomarkers it stays fixed at 1,
the identification constraint for a probit latent variable); only `D` is
re-estimated.

## Usage

``` r
longitudinalSubVarCopula(
  thetaLong,
  l,
  biomarker_type,
  thresholds,
  tol.em = 1e-04,
  max.iter = 100,
  verbose = FALSE
)
```
