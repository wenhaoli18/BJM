# Restricted cubic spline basis (Royston–Parmar)

Columns \\x\\ and, for each internal knot \\k_j\\, \\(x - k_j)\_+^3 -
\lambda_j (x - k\_{min})\_+^3 - (1 - \lambda_j) (x - k\_{max})\_+^3\\,
\\\lambda_j = (k\_{max} - k_j) / (k\_{max} - k\_{min})\\: a cubic spline
that is linear beyond the boundary knots. With `deriv = TRUE`, the
derivative of each column in \\x\\.

## Usage

``` r
rp_basis(x, knots, deriv = FALSE)
```

## Arguments

- x:

  Values (log times).

- knots:

  Boundary and internal knots, sorted, boundary first and last.

- deriv:

  Return the derivative instead.
