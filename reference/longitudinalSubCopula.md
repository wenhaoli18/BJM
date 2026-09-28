# Fit a multivariate longitudinal sub-model with mixed continuous/ordinal biomarkers (Gaussian-copula extension)

Internal counterpart to
[`longitudinalSubGaussian()`](https://liwh0904.github.io/BJM/reference/longitudinalSubGaussian.md),
used by the
[`longitudinalSub`](https://liwh0904.github.io/BJM/reference/longitudinalSub.md)
dispatcher whenever at least one biomarker is categorical/ordinal. Each
biomarker is still fit on its own first, exactly following
[`longitudinalSubGaussian()`](https://liwh0904.github.io/BJM/reference/longitudinalSubGaussian.md)'s
IFM strategy: continuous biomarkers via
[`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) (byte-identical
code path to
[`longitudinalSubGaussian()`](https://liwh0904.github.io/BJM/reference/longitudinalSubGaussian.md)),
ordinal biomarkers via
[`fit_marginal_ordinal()`](https://liwh0904.github.io/BJM/reference/fit_marginal_ordinal.md)
([`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html) with a
probit link). The resulting per-biomarker fixed effects, thresholds, and
random-effects covariance blocks seed the same `D`/`beta`/`sigma2`
bookkeeping
[`longitudinalSubGaussian()`](https://liwh0904.github.io/BJM/reference/longitudinalSubGaussian.md)
builds, which is then handed to
[`longitudinalSubVarCopula()`](https://liwh0904.github.io/BJM/reference/longitudinalSubVarCopula.md)
(rather than
[`longitudinalSubVar()`](https://liwh0904.github.io/BJM/reference/longitudinalSubVar.md))
to jointly re-estimate `D` across all biomarkers, including the ordinal
ones' latent Gaussian scores.

## Usage

``` r
longitudinalSubCopula(
  data_fit_all,
  long_sub_fixed,
  long_sub_random,
  biomarker_type,
  M
)
```
