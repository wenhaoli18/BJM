# Vertical line at a patient's observed event time

Helper for
[`riskPlot()`](https://wenhaoli18.github.io/BJM/reference/riskPlot.md):
a `geom_vline()` at the patient's survival time, coloured by event type
(black for 0, red for 1) with competing risks and red otherwise – or
`NULL` (nothing drawn) when that time is unknown, as for a new patient.
A missing time used to be passed to `geom_vline()` as is, and printing
the plot failed with "Discrete values supplied to continuous scale".

## Usage

``` r
event_time_vline(
  data_predict_all_pre,
  survival_variable,
  event_type_variable = NULL
)
```

## Arguments

- data_predict_all_pre:

  The patient's prediction data (list of data frames).

- survival_variable:

  Name of the survival-time variable.

- event_type_variable:

  Name of the event-type variable, or `NULL`.

## Value

A ggplot2 layer or `NULL`.
