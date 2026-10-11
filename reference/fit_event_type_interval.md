# Event-type model for an interval-censored outcome

The competing-risks part of
[`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md):
a logistic model for the event type among subjects whose event was
observed (finite `R`). The event type is taken as known once the event
is detected, but the event time is not, so a model that uses it needs
the `event_time` column filled for those subjects
([`fitIntervalBJM()`](https://wenhaoli18.github.io/BJM/reference/fitIntervalBJM.md)
fills it with imputed times and refits this model every iteration).

## Usage

``` r
fit_event_type_interval(
  data,
  form_marginal_surv,
  form_conditional_cr,
  event_time
)
```
