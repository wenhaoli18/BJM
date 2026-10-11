# Proportional hazards model for interval-censored event times

Internal fitter behind
[`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md)
when the outcome is interval-censored
(`Surv(L, R, type = "interval2")`).
[`coxph`](https://rdrr.io/pkg/survival/man/coxph.html) cannot fit such
data, so the marginal model is a proportional hazards model \\H(t \mid
x) = H_0(t) e^{x\beta}\\ with a parametric baseline, fit by maximum
likelihood. Each subject contributes \\S(L) - S(R)\\ for \\L \< T \le
R\\, \\S(L)\\ when right-censored (\\R = \infty\\), \\1 - S(R)\\ when
left-censored, and \\h(T) S(T)\\ for an exactly observed \\T\\. Two
baselines:

- `"spline"` (default):

  Royston–Parmar: \\\log H_0(t)\\ is a natural cubic spline in \\\log
  t\\ with `df` degrees of freedom (`df - 1` internal knots at quantiles
  of the log interval endpoints, boundary knots at their extremes).
  `df = 1` is a Weibull model.

- `"piecewise"`:

  a piecewise-constant hazard on `df` pieces cut at quantiles of the
  interval endpoints.

The spline is the default because interval-censored data say little
about the shape of the hazard *within* a visit interval, and that shape
is exactly what
[`imputeEventTime()`](https://wenhaoli18.github.io/BJM/reference/imputeEventTime.md)
draws \\T\\ from: the piecewise-constant hazard is flat there, which in
simulations biased the longitudinal sub-model fit by
[`fitIntervalBJM()`](https://wenhaoli18.github.io/BJM/reference/fitIntervalBJM.md)
when the true hazard was increasing, while the smooth spline carries the
trend across neighbouring intervals into each one.

The baseline cumulative hazard is then tabulated on a fine grid in the
same `hazard`/`time` layout as
[`basehaz`](https://rdrr.io/pkg/survival/man/basehaz.html), so the
prediction code can use it exactly as it uses the Breslow estimate of a
right-censored fit (see `survival_accessors`).

## Usage

``` r
icphFit(
  formula,
  data,
  event_time,
  baseline = c("spline", "piecewise"),
  df = 3,
  lower_tail = c("spline", "weibull")
)
```

## Arguments

- formula:

  `Surv(L, R, type = "interval2") ~ covariates`. A
  [`strata()`](https://rdrr.io/pkg/survival/man/strata.html) term gives
  each stratum its own baseline (spline knots or cut points from that
  stratum's interval endpoints), with the regression coefficients
  shared, as in a stratified Cox model.

- data:

  One row per subject.

- event_time:

  Name of the event-time column the longitudinal sub-model uses.

- baseline:

  `"spline"` or `"piecewise"`.

- df:

  Spline degrees of freedom, or number of pieces.

- lower_tail:

  Spline only: how \\\log H_0\\ is extrapolated before the earliest
  interval endpoint (the lower boundary knot), where it never enters the
  likelihood. `"spline"` (default) continues the spline's own linear
  tail; `"weibull"` uses the slope of the nested Weibull fit (`df = 1`)
  instead. In simulations where a third to a half of the events fell
  before the first visit, `"weibull"` made
  [`fitIntervalBJM()`](https://wenhaoli18.github.io/BJM/reference/fitIntervalBJM.md)'s
  estimates less variable but biased them when the true early hazard was
  not Weibull-shaped (Gompertz), so it is kept only for sensitivity
  analysis: the shape there is not identified by interval-censored data,
  and every choice is an assumption.

## Value

An object of class `"icph.BJM"`.
