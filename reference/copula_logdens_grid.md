# A patient's copula log densities at every grid point

Shared by
[`conditionalYTCopula()`](https://wenhaoli18.github.io/BJM/reference/conditionalYTCopula.md),
[`conditionalYDTCopula()`](https://wenhaoli18.github.io/BJM/reference/conditionalYDTCopula.md),
[`conditionalYTBioCopula()`](https://wenhaoli18.github.io/BJM/reference/conditionalYTBioCopula.md)
and
[`conditionalYDTBioCopula()`](https://wenhaoli18.github.io/BJM/reference/conditionalYDTBioCopula.md).
The patient's covariance, observed values and threshold bounds do not
depend on the event time, so the design
([`build_conditional_design_copula()`](https://wenhaoli18.github.io/BJM/reference/build_conditional_design_copula.md))
is built once per candidate value of the predicted biomarker instead of
once per grid point, and the means at all grid points come from
[`copula_mean_grid()`](https://wenhaoli18.github.io/BJM/reference/copula_mean_grid.md).
All densities are evaluated at once.

Predicting a continuous biomarker with the Monte Carlo box probability
(two or more ordinal random effects, see
[`mixed_density_copula_prep()`](https://wenhaoli18.github.io/BJM/reference/mixed_density_copula_prep.md))
needs it at every grid point for every candidate value, which made this
the slowest prediction in BJM. The candidate value only shifts the
ordinal latent scores' conditional mean, so the log box probability is a
smooth function of it: it is computed at `n_exact` candidate values
spread over the grid, and a cubic spline gives it at the others, while
the continuous biomarkers' density is computed exactly at every value.
In checks the difference from computing every value was no larger than
the Monte Carlo error itself (around \\10^{-4}\\), at 4-8 times the
speed. With one ordinal random effect, the quadrature is fast and every
value is computed exactly.

## Usage

``` r
copula_logdens_grid(
  data_num_i_list,
  rep_num_i_list,
  long_fit_all,
  l_i,
  survival_variable,
  time_variable,
  survival_variable_all,
  survival_trans_function,
  sigma.longitudinal,
  event_type_variable = NULL,
  bio = NULL,
  n_exact = 15
)
```

## Arguments

- data_num_i_list, rep_num_i_list:

  The patient's rows per biomarker.

- sigma.longitudinal:

  Residual SDs (1 for an ordinal biomarker).

- event_type_variable:

  `NULL` without competing risks; otherwise the event-type column,
  evaluated at 1 and at 0.

- bio:

  NULL, or a list with `index`, `name`, `row` and `values`: the
  predicted biomarker, its response column, the row of the predicted
  measurement and its candidate values.

- n_exact:

  Number of candidate values computed exactly when interpolating (see
  Details).

## Value

A list with one element per event type (1, then 0; just one without
competing risks), each a `length(l_i)` x (number of candidate values, 1
without `bio`) matrix of log densities.
