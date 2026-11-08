# The "unlabeled" method of transportBJM()

Estimates the log baseline-hazard shift on the target histories and,
when `data_source` is given, on the source histories too, and returns
the difference (see
[`transportBJM()`](https://wenhaoli18.github.io/BJM/reference/transportBJM.md)).

## Usage

``` r
transport_unlabeled(
  source_fit,
  data_predict_all,
  long_fit_all,
  time_variable,
  survival_variable_all,
  survival_trans_function,
  landmark,
  n_grid,
  upper,
  alpha_range,
  data_source,
  cap_event_time
)
```
