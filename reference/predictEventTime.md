# Predictive distribution of the event time

**Experimental.** Each patient's whole predictive distribution of the
event time \\T\\ given the biomarker history and \\T \> s\\, \\p(t \mid
Y, T \> s) \propto f(Y \mid t) f(t)\\ for \\t \> s\\, in one pass: the
cumulative distribution function at `prediction_time + horizons` and its
quantiles. The CDF at a horizon is
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)'s
risk (summed over the event types under competing risks), but computing
the whole curve at once is much faster than calling
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
once per horizon. The quantiles give prediction intervals for \\T\\
itself, e.g. the 10\\ time the patient outlives with probability 0.9.

Beyond the last follow-up of the survival sub-model the distribution
rests on its extrapolated baseline hazard, so upper quantiles that fall
there are only as good as that tail (see
[`survivalSubExternal`](https://wenhaoli18.github.io/BJM/reference/survivalSubExternal.md)).
Within follow-up the Breslow baseline puts the probability at the
observed event times, so the CDF rises in steps there; quantiles are
interpolated within the grid cell where the CDF crosses the probability.

## Usage

``` r
predictEventTime(
  data_predict_all,
  long_fit_all,
  survival_fit_all,
  prediction_time,
  horizons,
  time_variable,
  survival_variable_all = NULL,
  survival_trans_function = NULL,
  probs = c(0.1, 0.5, 0.9),
  n_grid = 200,
  n_tail = 100
)
```

## Arguments

- data_predict_all, long_fit_all, survival_fit_all, time_variable,
  survival_variable_all, survival_trans_function:

  As in
  [`predictRisk`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
  (right-censored survival fits only).

- prediction_time:

  The landmark \\s\\.

- horizons:

  Times after `prediction_time` at which to evaluate the CDF.

- probs:

  Probabilities of the quantiles.

- n_grid:

  Grid cells between `prediction_time` and the last horizon or the
  survival sub-model's last follow-up, whichever is later; `n_tail` grid
  cells beyond it.

- n_tail:

  See `n_grid`.

## Value

An object of class `"predictEventTime.BJM"`, a list with `cdf` (patients
x horizons, \\P(T \le s + h \mid Y, T \> s)\\); with competing risks,
`cif`, a list of two such matrices, one per event type (0, 1; they sum
to `cdf`); `quantile` (patients x `probs`, absolute times);
`prediction_time`; `horizons`.

## Examples

``` r
# \donttest{
data(pbc3)
surv <- survivalSub(pbc3[!duplicated(pbc3$id), ], Surv(years, status3) ~ age + sex, NULL)
long <- longitudinalSub(pbc3[pbc3$status3 == 1, ], list(serBilir ~ year + age + years),
                        list(~ year | id))
history <- pbc3[pbc3$id %in% c(2, 7) & pbc3$year <= 3, ]
history$years <- NA
out <- predictEventTime(history, long, surv, prediction_time = 3, horizons = 1:5,
                        time_variable = "year")
out$cdf
#>           h1        h2        h3        h4        h5
#> 2 0.05116750 0.1195239 0.2082627 0.3702554 0.4548032
#> 7 0.03034362 0.0748631 0.1385939 0.2693125 0.3448145
out$quantile
#>        10%      50%      90%
#> 2 4.621278 8.441462 11.17829
#> 7 5.364772 8.990141 13.87635
# }
```
