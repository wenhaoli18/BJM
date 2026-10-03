# Set a patient's future covariates

For each row of `future_rows` (one per simulated time), the covariates
of the last row of `future_covariates` (for this patient, if it has an
id column) at or before that time; a covariate it does not set keeps its
value in `future_rows`.

## Usage

``` r
apply_future_covariates(
  future_rows,
  future_covariates,
  time_variable,
  id,
  patient_id
)
```
