# Slope of the extrapolated baseline cumulative hazard

The constant hazard
[`cumulative_baseline_at()`](https://wenhaoli18.github.io/BJM/reference/cumulative_baseline_at.md)
uses past the last tabulated time: the slope of the least-squares line
through the whole table, or 0 if that slope is negative or undefined (a
table with a single time), so the cumulative hazard never decreases.

## Usage

``` r
extrapolation_slope(cum_basehaz)
```

## Arguments

- cum_basehaz:

  A data frame with columns `hazard` and `time`, ordered by time.

## Value

A single non-negative number.
