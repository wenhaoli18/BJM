# Build the per-patient longitudinal design pieces for a mixed continuous/ordinal (Gaussian-copula) joint model

Copula-aware counterpart to
[`build_conditional_design()`](https://wenhaoli18.github.io/BJM/reference/build_conditional_design.md)
(see `R/conditionalDesign.R`), for use whenever
`long_fit_all$biomarker_type` contains at least one `"ordinal"`
biomarker (see
[`longitudinalSubCopula()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSubCopula.md)).
Continuous biomarkers' rows are handled identically to
[`build_conditional_design()`](https://wenhaoli18.github.io/BJM/reference/build_conditional_design.md):
their exact observed value is recorded. For an ordinal biomarker's rows,
the exact underlying latent Gaussian score is unobserved – only the pair
of cumulative-link thresholds (`long_fit_all$thresholds[[m]]`)
bracketing the observed category is known (the classical Albert & Chib
probit data-augmentation representation, same convention as
[`impute_latent_ordinal()`](https://wenhaoli18.github.io/BJM/reference/impute_latent_ordinal.md)
in `R/copulaLongitudinal.R`: latent score `Z = eta + epsilon`,
`epsilon ~ N(0, 1)`, category `k` iff
`full_alpha[k] < Z <= full_alpha[k+1]` with
`full_alpha = c(-Inf, alpha_m, Inf)`, thresholds used directly with no
sign flip).

The random-effects design (`A_i`) and `Sigma_all` construction is
unchanged in form from
[`build_conditional_design()`](https://wenhaoli18.github.io/BJM/reference/build_conditional_design.md)
– marker type only affects `sigma.longitudinal` (already fixed at `1`
for ordinal markers by the caller, the probit identification constraint;
see
[`longitudinalSubCopula()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSubCopula.md)),
which flows in exactly the same way for both continuous and ordinal
markers.

## Usage

``` r
build_conditional_design_copula(
  rep_num_i_list,
  data_num_i_list,
  lfit,
  Sigma,
  sigma.longitudinal,
  time_variable,
  n_longitudinal,
  biomarker_type,
  long_sub_fixed,
  thresholds
)
```

## Value

A list with `y_all_vec` (stacked observed values; `NA` at ordinal-marker
rows), `alpha_lower_vec`/`alpha_upper_vec` (stacked threshold bounds;
`NA` at continuous-marker rows), `row_marker_type` (stacked per-row
`"continuous"`/`"ordinal"` label), `parameter_matrix` (block-diagonal
fixed-effect coefficients, `lfit[[i]]$coefficients$fixed` for continuous
markers, `lfit[[i]]$beta` – no intercept – for ordinal markers), and
`Sigma_all`.
