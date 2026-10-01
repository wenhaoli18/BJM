# Time by which conditional survival falls below a threshold

Helper for
[`integration_upper_bound()`](https://wenhaoli18.github.io/BJM/reference/integration_upper_bound.md):
for each linear predictor in `lp`, the earliest time \\t \> s\\ with
\\\exp(-(H_0(t) - H_0(s)) e^{lp}) \<\\ `tail_prob`, where \\H_0\\ is the
tabulated baseline cumulative hazard and, past its last time, the same
least-squares line
[`marginalT()`](https://wenhaoli18.github.io/BJM/reference/marginalT.md)
extrapolates with. `Inf` if that line is not increasing.

## Usage

``` r
tail_time(cum_basehaz, lp, s, tail_prob)
```

## Arguments

- cum_basehaz:

  A data frame with columns `hazard` and `time`.

- lp:

  Linear predictors (`reference = "zero"`).

- s:

  The prediction time.

- tail_prob:

  The survival threshold.

## Value

A numeric vector, one time per element of `lp`.
