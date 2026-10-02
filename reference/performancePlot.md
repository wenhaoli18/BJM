# Plot predictive performance across landmark times

Evaluates how well the backward joint model's dynamic risk predictions
discriminate and are calibrated, at each of several landmark times. For
each landmark `s` in `prediction_time`, every subject in
`data_predict_all` still event-free at `s` gets a predicted risk from
[`predictRisk`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md),
using only their longitudinal history up to `s`, of an event in
`(s, s + horizon]`; this is compared with what actually happened, by

- AUC:

  the time-dependent (cumulative/dynamic) area under the ROC curve: the
  probability that a subject with an event in the window has a higher
  predicted risk than one still event-free at `s + horizon`. 0.5 is no
  better than chance, 1 is perfect.

- Brier score:

  the mean squared difference between the predicted risk and the
  observed event indicator in the window. Lower is better.

Subjects censored within the window are handled by inverse probability
of censoring weighting (IPCW), with the censoring distribution among the
subjects at risk at `s` estimated by Kaplan–Meier. With competing risks,
both measures are computed separately for each event type: the cases are
the subjects with that event type in the window, and for the Brier score
a subject with the other event type in the window counts as not having
the event.

Evaluating on the data the sub-models were fit on gives an optimistic
(apparent) performance; pass held-out validation data in
`data_predict_all` when available.

## Usage

``` r
performancePlot(
  data_predict_all,
  long_fit_all,
  survival_fit_all,
  prediction_time,
  horizon,
  time_variable,
  survival_variable_all,
  survival_trans_function,
  bandcount1 = "auto",
  bandcount2 = "auto"
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

A `ggplot` object with one panel per measure, the landmark time on the
horizontal axis. Its `data` element is a `data.frame` with one row per
landmark, measure (and event type), and columns `landmark`, `measure`,
`cause`, `value`, `n_at_risk` (subjects predicted) and `n_cases` (events
of that type observed in the window).

## Examples

``` r
# \donttest{
data(pbc3)
survival_fit_all <- survivalSub(pbc3[!duplicated(pbc3$id), ],
                                Surv(years, status3) ~ age + sex, NULL)
long_sub_fixed <- list("long1" = serBilir ~ year + age + sex + years)
long_sub_random <- list("long1" = ~ year | id)
long_fit_all <- longitudinalSub(list(pbc3[pbc3$status3 == 1, ]),
                                long_sub_fixed, long_sub_random)
trans <- survivalTrans(c(1, 3, 5, 7))

# apparent performance on the fitting data, 2-year window
performancePlot(pbc3, long_fit_all, survival_fit_all,
  prediction_time = c(1, 3, 5), horizon = 2, time_variable = "year",
  trans$survival_variable_all, trans$survival_trans_function)

# }
```
