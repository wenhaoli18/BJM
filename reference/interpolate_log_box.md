# Spline interpolation of a log box probability over candidate values

For
[`copula_logdens_grid()`](https://wenhaoli18.github.io/BJM/reference/copula_logdens_grid.md):
a cubic spline (end conditions from the end points' own cubic, not zero
curvature, which a log probability, curving like a parabola in its
tails, does not have) through the finite log probabilities at the
exactly computed values `x`, evaluated at `x_out`; outside the range of
the finite ones (where the probability underflowed to 0) the log
probability is `-Inf`.

## Usage

``` r
interpolate_log_box(x, y, x_out)
```

## Value

A vector, one value per element of `x_out`.
