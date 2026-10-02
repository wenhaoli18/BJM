# Restrict prediction data to patients still at risk

Shared helper for `predictRisk` and `dynamicPredictionBio`: drops rows
whose survival-time variable is below `prediction_time` from every
biomarker's data frame. Rows whose survival time is missing (e.g. a new
patient whose event time is not yet known) are kept: they are treated as
at risk. It is an error if that leaves no patient (it used to return an
empty result from
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
without comment, and fail with an unrelated error in
[`predictLongitudinal()`](https://wenhaoli18.github.io/BJM/reference/predictLongitudinal.md)).

## Usage

``` r
subset_at_risk(data_predict_all, survival_variable, prediction_time)
```

## Value

`data_predict_all`, filtered in place per element.
