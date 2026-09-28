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
survivalSub(data_survival_fitting, form_marginal_surv, form_conditional_cr)
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

## Value

An object of class `"survivalSub.BJM"`, a named list with elements:

- coxph_fit:

  The fitted [`coxph`](https://rdrr.io/pkg/survival/man/coxph.html)
  marginal survival model.

- form_marginal_surv:

  The `form_marginal_surv` formula, as supplied.

- glm_fit:

  The fitted [`glm`](https://rdrr.io/r/stats/glm.html) competing-risks
  (event type) model, or `NULL` if `form_conditional_cr` was not
  supplied.

- form_conditional_cr:

  The `form_conditional_cr` formula, as supplied (or `NULL`).

## Examples

``` r

data(pbc3)

data_survival_fitting =  pbc3[!duplicated(pbc3$id), ]

form_marginal_surv = Surv(years, status3) ~ age + sex
form_conditional_cr = NULL

survival_fit_all = survivalSub(data_survival_fitting, form_marginal_surv, 
                               form_conditional_cr)
```
