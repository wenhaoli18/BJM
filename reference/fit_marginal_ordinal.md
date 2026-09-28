# Fit the marginal model for a single ordinal/categorical biomarker

Fits a cumulative link mixed model (probit link) via
[`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html), used as
the marginal model for a categorical biomarker in the Gaussian-copula
extension of
[`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
(see
[`longitudinalSubCopula()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSubCopula.md)).
The probit link is used specifically (rather than the more common logit)
because it integrates directly with the latent-Gaussian-score
representation the copula extension relies on: the underlying continuous
latent variable implied by a probit cumulative link model is, by
construction, Gaussian, which is exactly the representation needed to
fold a categorical biomarker into the same random-effects covariance
structure as the continuous biomarkers.

This follows the classical two-stage IFM (Inference Functions for
Margins; Joe 2005) strategy already used throughout `BJM`: fit each
biomarker's marginal model separately first (here, exactly as
[`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) does for
continuous biomarkers), then re-estimate the joint random-effects
covariance across all biomarkers afterwards (see
[`longitudinalSubVarCopula()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSubVarCopula.md)).

## Usage

``` r
fit_marginal_ordinal(fixed_formula, random_formula, data)
```

## Value

A fitted [`ordinal::clmm`](https://rdrr.io/pkg/ordinal/man/clmm.html)
object.
