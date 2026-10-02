# Baseline cumulative hazard at a set of time points

Helper for `marginalT` and the integration-grid helpers: for each
element of `l_i`, the Breslow cumulative hazard as the right-continuous
step function it is – the value at the last tabulated time `<= l_i`, or
0 before the first. Past the table's last time it continues from the
last tabulated value with the slope of the least-squares line through
the whole table (a linear extrapolation of the cumulative hazard, i.e. a
constant hazard). It used to switch to that line itself, which does not
pass through the last tabulated value: the cumulative hazard jumped
there, and when it jumped down (a convex cumulative hazard, i.e. an
increasing hazard) the interval spanning the last time got a negative
probability mass. The line used to be tabulated on a fixed 0.005 grid
out to twice the integration upper limit and then looked up, which
assumed time was measured in years (with time in days the table ran to
millions of rows) and failed with "wrong sign in 'by'" when the upper
limit was below the last training time; it is now evaluated directly.
The table used to be read at the *nearest* tabulated time, which for a
point just before an event time took that event's jump early.

## Usage

``` r
cumulative_baseline_at(cum_basehaz, l_i)
```

## Arguments

- cum_basehaz:

  A data frame with columns `hazard` and `time`.

- l_i:

  Time points.

## Value

A numeric vector, one cumulative hazard per element of `l_i`.
