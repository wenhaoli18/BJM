# Interval endpoints from an interval-type `Surv` object

\\L \< T \le R\\ for every subject: \\R = \infty\\ when right-censored,
\\L = 0\\ when left-censored, \\L = R\\ when \\T\\ is observed exactly.

## Usage

``` r
interval_bounds(y)
```

## Arguments

- y:

  A `Surv` object of type `"interval"`.

## Value

A list with numeric vectors `L` and `R`.
