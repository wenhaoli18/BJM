# Assert that an object is a list of data.frame objects of the expected length

Shared input-validation helper for the `data_fit_all`/
`data_predict_all` arguments, which must be a list with one data.frame
per longitudinal outcome. When `allow_bare_df = TRUE`, a single bare
data.frame is also accepted, matching the auto-repeat convenience some
functions apply for a single data.frame shared across every outcome.

## Usage

``` r
assert_data_list(x, arg_name, n_expected, allow_bare_df = FALSE)
```
