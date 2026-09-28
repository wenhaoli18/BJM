# Restrict prediction data to patients still at risk

Shared helper for `predictRisk` and `dynamicPredictionBio`: drops rows
whose survival-time variable is below `prediction_time` from every
biomarker's data frame.

## Usage

``` r
subset_at_risk(data_predict_all, survival_variable, prediction_time)
```

## Value

`data_predict_all`, filtered in place per element.
