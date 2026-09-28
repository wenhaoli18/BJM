# Barnard & Rubin (1999) small-sample pooled degrees of freedom

Internal helper for
[`rubin_pool_scalar`](https://liwh0904.github.io/BJM/reference/rubin_pool_scalar.md).
Interpolates between the classical Rubin (1987) degrees of freedom
(`df_old`, which ignores the complete-data sample size and can therefore
exceed it) and the complete-data degrees of freedom (`dfcom`), weighted
by how much of the total variance is due to missingness (`lambda`).

## Usage

``` r
barnard_rubin_df(m, lambda, dfcom)
```
