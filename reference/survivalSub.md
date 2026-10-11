# Fit the survival sub-model

Fits the survival building block used by
[`predictRisk`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
and
[`dynamicPredictionBio`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md):
a marginal Cox proportional-hazards model (via
[`coxph`](https://rdrr.io/pkg/survival/man/coxph.html), with
`x = TRUE, y = TRUE` so the fit is self-contained for later prediction)
for the time-to-event outcome given `form_marginal_surv`, plus,
optionally, a logistic event-type/competing-risks model (via
[`glm`](https://rdrr.io/r/stats/glm.html) with `family = binomial`)
given `form_conditional_cr`. The competing-risks model is fit only among
subjects who experienced an event (i.e. whose censoring indicator is
non-zero), predicting which type of event occurred conditional on an
event having occurred. Together with
[`longitudinalSub`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md),
the object returned here forms the pair of sub-models that dynamic
prediction is built on.

## Usage

``` r
survivalSub(
  data_survival_fitting,
  form_marginal_surv,
  form_conditional_cr,
  event_time = NULL,
  baseline = c("spline", "piecewise"),
  df = NULL
)
```

## Arguments

- data_survival_fitting:

  One row per subject, containing the time-to-event outcome,
  censoring/event-type indicator, and any baseline covariates referenced
  in `form_marginal_surv` or `form_conditional_cr`.

- form_marginal_surv:

  A survival formula, e.g. `Surv(time, status) ~ covariates`, passed to
  [`coxph`](https://rdrr.io/pkg/survival/man/coxph.html) to fit the
  marginal event-time model.

- form_conditional_cr:

  An optional formula for the competing-risks (event-type) model, e.g.
  `event_type ~ covariates`, passed to
  [`glm`](https://rdrr.io/r/stats/glm.html) with `family = binomial`.
  Set to `NULL` when there is only a single event type (no competing
  risks).

- event_time:

  Only for an interval-censored outcome: the name of the event-time
  column that the longitudinal sub-model's formulas use (the exact time
  is unknown, so it is not part of `form_marginal_surv`).

- baseline:

  Only for an interval-censored outcome: the baseline hazard, `"spline"`
  (default; Royston–Parmar spline for the log cumulative hazard in log
  time) or `"piecewise"` (piecewise-constant hazard).

- df:

  Only for an interval-censored outcome: spline degrees of freedom
  (default 3) or number of pieces (default 6).

## Value

An object of class `"survivalSub.BJM"`, a named list with elements:

- coxph_fit:

  The fitted [`coxph`](https://rdrr.io/pkg/survival/man/coxph.html)
  marginal survival model (`ic_fit`, an interval-censored PH fit, in its
  place for an interval-censored outcome).

- form_marginal_surv:

  The `form_marginal_surv` formula, as supplied.

- glm_fit:

  The fitted [`glm`](https://rdrr.io/r/stats/glm.html) competing-risks
  (event type) model, or `NULL` if `form_conditional_cr` was not
  supplied.

- form_conditional_cr:

  The `form_conditional_cr` formula, as supplied (or `NULL`).

## Details

**Interval censoring (experimental).** When `form_marginal_surv` is
`Surv(L, R, type = "interval2") ~ covariates` – the event is only known
to lie in `(L, R]`, with `R = NA` or `Inf` for right-censored subjects –
the marginal model is a proportional hazards model with a smooth spline
(or piecewise-constant) baseline hazard fit by maximum likelihood (see
[`icphFit()`](https://wenhaoli18.github.io/BJM/reference/icphFit.md)),
stored as `ic_fit` in place of `coxph_fit`. A
[`strata()`](https://rdrr.io/pkg/survival/man/strata.html) term gives
each stratum its own baseline hazard. With competing risks, the event
type (0/1) is taken as known once the event is detected, and
`form_conditional_cr` is fit among subjects with an observed event
(finite `R`); if it uses the event time, that column must be filled for
them, as
[`fitIntervalBJM`](https://wenhaoli18.github.io/BJM/reference/fitIntervalBJM.md)
does with imputed times. Ordinary right-censored `Surv(time, status)`
outcomes are fit with `coxph` exactly as before.

Events that happened before a subject's first visit (`L = 0`) only say
that the event fell between time 0 and that visit. Before the earliest
visit in the data there are no interval endpoints, so the baseline
hazard there is extrapolated, not estimated, and
[`fitIntervalBJM`](https://wenhaoli18.github.io/BJM/reference/fitIntervalBJM.md)
places these subjects' event times according to that extrapolation.
[`print()`](https://rdrr.io/r/base/print.html) reports how many events
are of this kind: when they are a large share of all events (in
simulations, around a third or more), results can depend noticeably on
the assumed shape of the early hazard.

## Examples

``` r

data(pbc3)

data_survival_fitting =  pbc3[!duplicated(pbc3$id), ]

form_marginal_surv = Surv(years, status3) ~ age + sex
form_conditional_cr = NULL

survival_fit_all = survivalSub(data_survival_fitting, form_marginal_surv, 
                               form_conditional_cr)
```
