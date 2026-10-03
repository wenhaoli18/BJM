# Fixed-effects means at every candidate event time

For each biomarker's rows, \\X(l)\beta\\ with the survival-time columns
(and the event type, with competing risks) set to each element of
`l_grid`, stacked across biomarkers.

## Usage

``` r
stacked_mean_on_grid(
  rows,
  long_fit_all,
  l_grid,
  d,
  survival_variable,
  event_type_variable,
  survival_variable_all,
  survival_trans_function
)
```

## Value

A matrix with one row per stacked observation and one column per element
of `l_grid`.
