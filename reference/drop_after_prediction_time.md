# Drop longitudinal measurements taken after the prediction time

Shared helper for `predictRisk`, `dynamicPredictionBio`, and
`dynamicPredictionBioAll`: a dynamic prediction at `prediction_time` may
only condition on the history observed up to `prediction_time`, so rows
of `data_predict_all` whose `time_variable` is later than that are
removed, with a warning saying how many. (`predictPlot`/`riskPlot`
already truncate per landmark time before calling these, so they never
trigger the warning.) Rows with a missing `time_variable` are kept, as
before.

## Usage

``` r
drop_after_prediction_time(data_predict_all, time_variable, prediction_time)
```

## Value

`data_predict_all`, filtered per element.
