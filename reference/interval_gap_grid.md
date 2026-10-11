# Integration grid between the last negative visit and the prediction time

Interval-censored fits only. With interval censoring an event is only
detected at the next visit, so at `prediction_time` \\s\\ a patient is
known to be event-free only up to their last visit \\V \le s\\ (taken to
be the last time in their biomarker history up to \\s\\), not up to
\\s\\ itself. Conditioning on \\T \> V\\ instead of \\T \> s\\ adds the
stretch \\(V, s\]\\ to the denominator of every prediction. This builds,
for each patient, `n` equal intervals tiling their own \\(V, s\]\\ – the
stretch differs per patient, so the grid does too – and the survival
model's probability of each.

Returns `NULL` for a right-censored fit (an event is observed when it
happens, so \\V = s\\) and when every patient was seen at \\s\\, so
callers keep their right-censored code path unchanged.

## Usage

``` r
interval_gap_grid(
  data_predict_all,
  long_fit_all,
  survival_fit_all,
  prediction_time,
  time_variable,
  n
)
```

## Arguments

- data_predict_all:

  At-risk prediction data (list of data frames).

- long_fit_all:

  Output of
  [`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md).

- survival_fit_all:

  Output of
  [`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md).

- prediction_time:

  The prediction time \\s\\.

- time_variable:

  Visit-time column.

- n:

  Grid intervals per patient.

## Value

`NULL`, or a list with `last_visit` (\\V\\, one per patient), `mids` and
`S` (`n` x patients matrices of interval midpoints and probabilities;
zero probability for a patient with \\V = s\\), and `eval(fun)`, which
calls `fun(data_j, l_i)` once per patient with a gap (`data_j` the
patient's rows, `l_i` their midpoints) and binds the resulting
one-column results into an `n` x patients matrix (filled with `fill`,
`-Inf` by default for log densities, for patients without a gap). `fun`
may also return a (nested) list of such columns, in which case
[`eval()`](https://rdrr.io/r/base/eval.html) returns the matching list
of matrices.
