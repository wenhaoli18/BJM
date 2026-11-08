# Warn once when a prediction reaches beyond the survival model's follow-up

Called by
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
and
[`predictEventTime()`](https://wenhaoli18.github.io/BJM/reference/predictEventTime.md):
past the last follow-up time of the survival sub-model its baseline
hazard is extrapolated, and risks there can be far off (see
[`fitImputedBJM`](https://wenhaoli18.github.io/BJM/reference/fitImputedBJM.md)).
Warns once per session; set `options(BJM.check_extrapolation = FALSE)`
to silence it.

## Usage

``` r
check_beyond_followup(survival_fit_all, end)
```

## Arguments

- end:

  The last time predicted (`prediction_time + horizon`).
