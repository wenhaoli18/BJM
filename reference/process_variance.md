# Factor one patient's joint longitudinal covariance for conditionalYTBio()/conditionalYDTBio()

Factor one patient's joint longitudinal covariance for
conditionalYTBio()/conditionalYDTBio()

## Usage

``` r
process_variance(
  num_i,
  time_new,
  bio_i,
  data_predict_all,
  long_fit_all,
  time_variable
)
```

## Value

The result of
[`cov_factor()`](https://wenhaoli18.github.io/BJM/reference/cov_factor.md),
or `NA` if the patient has no rows for some biomarker.
