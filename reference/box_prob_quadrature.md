# Ordinal box probabilities by adaptive Gauss-Hermite quadrature

For each column \\k\\, \\\log P(lower \< m + l u + e \le upper)\\ with
\\u \sim N(0, 1)\\ and \\e \sim N(0, I)\\ independent: \\\log \int
\prod_j P_j(u) \phi(u) du\\, where \\P_j\\ is a univariate normal
probability. The integrand is log-concave in \\u\\; Newton's method
finds its mode and the curvature there, which place an `n_nodes`-point
Gauss-Hermite rule, all columns at once. With 20 nodes the relative
error was below \\10^{-5}\\ in every check.

## Usage

``` r
box_prob_quadrature(m, lower, upper, l, n_nodes = 20)
```

## Arguments

- m, lower, upper:

  Matrices, one row per ordinal observation and one column per
  evaluation.

- l:

  Output of
  [`latent_factor()`](https://wenhaoli18.github.io/BJM/reference/latent_factor.md).

- n_nodes:

  Number of quadrature nodes.

## Value

A vector of log probabilities.
