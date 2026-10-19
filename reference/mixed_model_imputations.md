# Multiple imputation from the multivariate linear mixed model

`method = "mixed"` in
[`imputeLongitudinal`](https://wenhaoli18.github.io/BJM/reference/imputeLongitudinal.md).
Fits the model
[`longitudinalSub`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
fits – each biomarker's linear mixed model, with the random effects of
all biomarkers jointly normal with covariance `Sigma_fit` – and draws
each subject's missing values from their normal distribution given
everything observed for that subject: all of its visits and all
biomarkers. For each completion the model's parameters are first drawn
from an approximation to their posterior, so that the completions also
reflect the uncertainty in the fit (proper multiple imputation): the
fixed effects from their estimated sampling distribution, each residual
variance from a scaled inverse chi-squared distribution, and `Sigma_fit`
from an inverse Wishart distribution with as many degrees of freedom as
subjects.

## Usage

``` r
mixed_model_imputations(
  aligned,
  data_fit_all,
  long_sub_fixed,
  long_sub_random,
  id,
  resp_mat,
  mask,
  n_imputations
)
```

## Value

A list of `n_imputations` matrices like `resp_mat`, with the missing
cells filled.
