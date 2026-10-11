# Baselines for `icphFit()`

Each returns a list with `n_par`, `cumhaz(t, g)` (baseline cumulative
hazard, 0 at `t <= 0`, `Inf` at `Inf`), `loghaz(t, g)` (log baseline
hazard, for exactly observed times), `monotone(g)` (whether `cumhaz` is
increasing), `start(rate0)` (parameters of a constant hazard `rate0`)
and `knot_times`.

## Usage

``` r
rp_baseline(endpoints, df, tail_slope = NULL)

piecewise_baseline(endpoints, df)
```

## Arguments

- endpoints:

  Finite, positive interval endpoints.

- df:

  Spline degrees of freedom / number of pieces.

## Functions

- `rp_baseline()`: Royston–Parmar spline for \\\log H_0\\ in \\\log t\\.
  Parameters: intercept, then one per basis column; the basis columns
  are scaled by their standard deviation over `endpoints` so the
  optimizer sees comparable scales.

- `piecewise_baseline()`: Piecewise-constant hazard; parameters are the
  log hazards of the pieces.
