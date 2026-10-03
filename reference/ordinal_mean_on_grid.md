# Latent-score means of an ordinal biomarker at every candidate event time

\\X(l)\beta\\ for an ordinal biomarker's rows, with the design built the
way
[`build_LME_indi_matrix_copula()`](https://wenhaoli18.github.io/BJM/reference/build_LME_indi_matrix_copula.md)
builds it (the columns of `lfit[[i]]$beta`, which has no intercept: the
thresholds play that role).

## Usage

``` r
ordinal_mean_on_grid(
  data_i,
  i,
  long_fit_all,
  xlev_i,
  l_unique,
  survival_variable,
  survival_variable_all,
  survival_trans_function
)
```

## Value

A matrix with one row per row of `data_i` and one column per element of
`l_unique`.
