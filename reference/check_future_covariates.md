# Validate `future_covariates`

Checks that it is a data frame with a non-missing time column, that its
id column (if any) names patients being simulated, and that it does not
set an outcome. Returns it ordered by time.

## Usage

``` r
check_future_covariates(
  future_covariates,
  time_variable,
  id,
  patient_ids,
  outcome_variables
)
```
