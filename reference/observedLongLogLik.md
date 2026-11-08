# Observed-data log-likelihood of the longitudinal sub-model

**Experimental.** Compares longitudinal sub-models fitted to the same
right-censored data, e.g. different forms of the event-time effect in
`long_sub_fixed`, or a fit on the subjects with an event against
[`fitImputedBJM`](https://wenhaoli18.github.io/BJM/reference/fitImputedBJM.md).
A subject with an event contributes \\\log f(Y \mid T, D)\\ at the
observed \\T\\ (and event type); a censored subject contributes \\\log
\int_C^\infty \sum_d f(Y \mid t, d) P(d \mid t) \\ dF(t) - \log S(C)\\,
so the censored subjects count too. The factors \\f(T)\\ and \\P(D \mid
T)\\ of the events do not involve the longitudinal sub-model and are
left out, so only differences between fits sharing `survival_fit_all`
are meaningful. Higher is better; the fits compared should have the same
number of parameters, or the difference be judged with that in mind.

## Usage

``` r
observedLongLogLik(
  fit,
  data_survival,
  data_fit_all,
  time_variable,
  survival_variable_all = NULL,
  survival_trans_function = NULL,
  n_grid = 40
)
```

## Arguments

- fit:

  A `"fitImputedBJM"` object (its log-likelihood is averaged over the
  imputations' fits), or a list with `long_fit_all` and
  `survival_fit_all`.

- data_survival, data_fit_all:

  The data the fit was made on, as in
  [`fitImputedBJM`](https://wenhaoli18.github.io/BJM/reference/fitImputedBJM.md):
  every subject, with the observed event or censoring time.

- time_variable, survival_variable_all, survival_trans_function:

  As in
  [`predictRisk`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md).

- n_grid:

  Grid points per censored subject.

## Value

The log-likelihood (a number).
