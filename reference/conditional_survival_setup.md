# Per-patient pieces of the marginal survival model

Shared by
[`integration_upper_bound()`](https://wenhaoli18.github.io/BJM/reference/integration_upper_bound.md)
and
[`prepare_infinity_grid()`](https://wenhaoli18.github.io/BJM/reference/prepare_infinity_grid.md):
each at-risk patient's Cox linear predictor (`reference = "zero"`, `NA`
if a covariate is missing), and the patients grouped by the baseline
cumulative hazard they use (one group, or one per stratum of a
stratified model).

## Usage

``` r
conditional_survival_setup(data_predict_all, long_fit_all, survival_fit_all)
```

## Value

A list with `lp`, `groups` (each a list with `cum_basehaz` and
`patients`, indices into `lp`), and `last_time` (the last time in the
data
[`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md)
was fit on).
