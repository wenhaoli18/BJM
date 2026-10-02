# Predicted risks and observed outcomes at each landmark time

Shared by
[`performancePlot()`](https://wenhaoli18.github.io/BJM/reference/performancePlot.md)
and
[`calibrationPlot()`](https://wenhaoli18.github.io/BJM/reference/calibrationPlot.md):
validates their common arguments and, for each landmark `s`, calls
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
on the subjects event-free at `s`, with their history up to `s` and
their outcome hidden.

## Usage

``` r
landmark_predictions(
  data_predict_all,
  long_fit_all,
  survival_fit_all,
  prediction_time,
  horizon,
  time_variable,
  survival_variable_all,
  survival_trans_function,
  bandcount1,
  bandcount2
)
```

## Arguments

- data_predict_all:

  The evaluation data, in the same format as for
  [`predictRisk`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md):
  a list of long-format `data.frame`s, one per longitudinal outcome (or
  a single `data.frame` used for all). Unlike for
  [`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md),
  it must also contain each subject's observed outcome: the survival
  time and status variables of `survival_fit_all`'s Cox formula and,
  with competing risks, the event type variable of its
  `form_conditional_cr`. The outcome is taken from the first data frame
  and is hidden from
  [`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md).

- long_fit_all:

  A `longitudinalSub.BJM` object from
  [`longitudinalSub`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md).

- survival_fit_all:

  A `survivalSub.BJM` object from
  [`survivalSub`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md).

- prediction_time:

  A vector of landmark times.

- horizon:

  The prediction horizon (a single positive number): the window after
  each landmark in which events are counted.

- time_variable:

  The name of the time variable in the linear mixed models.

- survival_variable_all, survival_trans_function:

  As for
  [`predictRisk`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md).

- bandcount1, bandcount2:

  As for
  [`predictRisk`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md);
  passed to it at each landmark.

## Value

A list with one element per landmark that had subjects to predict, each
a list with `landmark`, `outcome` (see
[`performance_outcome()`](https://wenhaoli18.github.io/BJM/reference/performance_outcome.md),
rows matching the predictions) and `risks` (a list with one vector of
predicted risks per event type).
