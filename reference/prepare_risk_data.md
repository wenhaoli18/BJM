# Prediction data for `predictRisk()`, keeping patients without measurements

Without any biomarker measurement up to `prediction_time`, the
longitudinal factor \\f(Y \mid T)\\ of the backward joint model is 1,
and the predicted risk is that of the survival sub-model alone given the
patient's covariates.
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
predicts such patients that way instead of refusing them. This runs the
usual
[`drop_after_prediction_time()`](https://wenhaoli18.github.io/BJM/reference/drop_after_prediction_time.md),
[`align_ordinal_levels()`](https://wenhaoli18.github.io/BJM/reference/align_ordinal_levels.md)
and
[`drop_missing_longitudinal()`](https://wenhaoli18.github.io/BJM/reference/drop_missing_longitudinal.md),
then puts back one row per patient left without any measurement: their
earliest row, with the biomarker set to `NA` and `time_variable` set to
their last visit up to `prediction_time` (`min(0, prediction_time)` if
there is none), which an interval-censored fit takes as the last time
they were known to be event-free. The covariates come from that row: the
survival sub-model's covariates are baseline values, one row per
patient. A warning names these patients.
[`split_no_measurement()`](https://wenhaoli18.github.io/BJM/reference/split_no_measurement.md)
separates them again; running this on its own output returns the same
data, as `"auto"` bandcount tuning does.

A patient with measurements of some biomarkers but not of others is
still reported as unpredictable, as before.

## Usage

``` r
prepare_risk_data(
  data_predict_all,
  long_fit_all,
  survival_fit_all,
  prediction_time,
  time_variable,
  survival_variable_all
)
```

## Value

`data_predict_all`, filtered, with the added rows in the position of
each patient's earliest row.
