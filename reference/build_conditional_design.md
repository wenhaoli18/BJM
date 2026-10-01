# Build the per-patient longitudinal design matrices

Shared helper for `conditionalYT` and `conditionalYDT`: constructs the
longitudinal outcome matrix, fixed-effect parameter matrix, and
random-effect variance-covariance pieces used by the conditional density
quadratic form, for a single patient.

## Usage

``` r
build_conditional_design(
  rep_num_i_list,
  data_num_i_list,
  lfit,
  Sigma,
  sigma.longitudinal,
  time_variable,
  n_longitudinal,
  long_sub_random
)
```

## Value

A list with `longitudinal_all_matrix`, `parameter_matrix`, `Sigma_all`,
`log_det_Var_cov_estep`, `Sigma_all_solve`, and `long_sigma_long`.
