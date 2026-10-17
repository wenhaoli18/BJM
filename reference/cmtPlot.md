# Plot conditional mean trajectories (CMT)

This function generates the Conditional Mean Trajectories (CMT) plot:
the mean biomarker trajectory of the patients whose event happened in
the window `(condi_time2event - event_window, condi_time2event]`. The
time from the first measurement up to `condi_time2event` is cut into
consecutive intervals of length `interval_time`; each point is the mean
of those patients' biomarker values measured within one interval,
plotted at the interval's end, and a straight line is fitted through the
points. With competing risks there is one trajectory per event type.

All times – `time_variable`, `survival_variable`, `condi_time2event`,
`event_window` and `interval_time` – are in the units of the data
(years, months, days, ...), and the time axis shows `time_variable` on
that scale, labelled with its name.

Only patients with an observed event enter: with competing risks, those
whose `event_type_variable` equals one of `event_type`; without, those
whose `status_variable` is not 0. Without competing risks and with
`status_variable = NULL`, patients are selected by their survival time
alone, so censored ones are included too, with a warning. An interval
without any measurement is left out rather than filled in.

## Usage

``` r
cmtPlot(
  data_plot_all,
  condi_time2event,
  event_type_variable,
  event_type,
  bio_variable,
  time_variable,
  survival_variable,
  interval_time = 1/12,
  id_variable = "id",
  status_variable = NULL,
  event_window = 1
)
```

## Arguments

- data_plot_all:

  A `data.frame` that includes the biomarker used for plotting. It is
  utilized to plot conditional mean trajectories (CMT).

- condi_time2event:

  Event time the patients are selected around. If `NULL`, it defaults to
  the midpoint of the observed range of `time_variable` in
  `data_plot_all`.

- event_type_variable:

  Competing risks variable indicator name. Set to NULL if there are no
  competing risks.

- event_type:

  A vector containing the values of `event_type_variable` for the event
  types to plot (at least two).

- bio_variable:

  Name of the biomarker variable used for plotting.

- time_variable:

  The name of time variable in linear mixed model.

- survival_variable:

  Name of the time-to-event outcomes variable.

- interval_time:

  Length of the intervals the follow-up is cut into, in the units of
  `time_variable`.

- id_variable:

  Name of the patient ID column in `data_plot_all`. Default is `"id"`.

- status_variable:

  Without competing risks, the name of the event indicator (0 =
  censored), so that only patients with an event are used. Ignored with
  competing risks.

- event_window:

  Width of the window of event times before `condi_time2event` that a
  patient's event must fall in, in the units of `survival_variable`.
  Default 1 (one year for data in years; use e.g. 12 for months or 365
  for days).

## Value

Conditional mean trajectories plot.

## Examples

``` r

# example without competing risks

data(pbc3)

pbc.cmt <- cmtPlot(data_plot_all = pbc3, condi_time2event = 5,
   event_type_variable = NULL, event_type = NULL,
   bio_variable = "serBilir", time_variable = "year",
   survival_variable = "years",
   interval_time = 1/12, status_variable = "status3"
)

pbc.cmt



# example with competing risks

data(pbc3)

pbc.cmt.cr <- cmtPlot(pbc3, condi_time2event = 5,
   event_type_variable = 'status4', event_type = c("0", "1"),
   bio_variable = "albumin", time_variable = "year",
   survival_variable = "years",
   interval_time = 1/4
)

pbc.cmt.cr

```
