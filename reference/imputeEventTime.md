# Draw interval-censored event times given the biomarker history

**Experimental.** One draw of each subject's event time from \\f(T \mid
L \< T \le R, Y, X) \propto f(Y \mid T, X) f(T \mid X)\\: \\f(Y \mid
T)\\ from the longitudinal sub-model on a grid of `n_grid` intervals
tiling \\(L, R\]\\
([`conditionalYT()`](https://wenhaoli18.github.io/BJM/reference/conditionalYT.md)),
times the survival model's probability of each grid interval
([`marginalT()`](https://wenhaoli18.github.io/BJM/reference/marginalT.md));
an interval is drawn with those weights and \\T\\ uniformly within it. A
right-censored subject (\\R = \infty\\) uses the same equal-probability
grid from \\L\\ out to the survival model's tail as
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md).
A subject with \\L = R\\ keeps \\T = L\\. Used by
[`fitIntervalBJM`](https://wenhaoli18.github.io/BJM/reference/fitIntervalBJM.md).

## Usage

``` r
imputeEventTime(
  data_fit_all,
  long_fit_all,
  survival_fit_all,
  bounds,
  time_variable,
  survival_variable_all = NULL,
  survival_trans_function = NULL,
  n_grid = 40,
  n_draws = 1
)
```

## Arguments

- data_fit_all:

  Longitudinal data, one data frame per biomarker.

- long_fit_all:

  Output of
  [`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md).

- survival_fit_all:

  Interval-censored output of
  [`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md).

- bounds:

  Data frame with columns `id`, `L`, `R` (and, with competing risks, the
  observed event type `D`).

- time_variable, survival_variable_all, survival_trans_function:

  As in
  [`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md).

- n_grid:

  Grid intervals per subject.

- n_draws:

  Independent draws per subject (the weights are computed once and
  reused).

## Value

A named numeric vector of draws, one per row of `bounds`; with
`n_draws > 1`, a matrix with one column per draw.
