# Gauss-Hermite quadrature rule for the standard normal

Nodes and weights (summing to 1) of the `k`-point rule for \\E f(X)\\,
\\X \sim N(0, 1)\\, by the Golub-Welsch eigenvalue method.

## Usage

``` r
gauss_hermite_normal(k)
```

## Arguments

- k:

  Number of nodes.

## Value

A list with `x` (nodes) and `w` (weights).
