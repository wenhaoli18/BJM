# Vertical line at the end of a patient's follow-up

Helper for
[`riskPlot()`](https://wenhaoli18.github.io/BJM/reference/riskPlot.md):
a vertical line at the patient's survival time, labelled in a legend –
solid, coloured by event type (as the risk curves: black for the first,
red for the second), for an observed event, and dashed grey for a
censored follow-up. Nothing is drawn when the time or the status is
unknown, as for a new patient or an interval-censored fit. The line used
to be drawn the same way (red, with no legend) for a censored patient as
for an event.

## Usage

``` r
event_time_vline(
  data_predict_all_pre,
  survival_fit_all,
  event_type_variable = NULL
)
```

## Arguments

- data_predict_all_pre:

  The patient's prediction data (list of data frames).

- survival_fit_all:

  The survival sub-model, whose
  [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) outcome gives
  the time and the status.

- event_type_variable:

  Name of the event-type variable, or `NULL`.

## Value

A list of ggplot2 components (possibly empty).
