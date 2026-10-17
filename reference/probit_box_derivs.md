# Log probability of an ordinal category and its derivatives

\\\log P(lo \< \eta + \epsilon \le hi)\\ for \\\epsilon \sim N(0, 1)\\,
with its first and second derivatives in \\\eta\\; evaluated on the side
of the mean where the interval lies, so a category far in a tail keeps
its precision.

## Usage

``` r
probit_box_derivs(eta, lo, hi)
```

## Arguments

- eta:

  Linear predictors; `lo`, `hi` the category bounds.

- lo, hi:

  Lower and upper thresholds (`-Inf`/`Inf` at the ends).

## Value

A list with `logP`, `g` (first derivative) and `h` (second derivative,
\\\le 0\\).
