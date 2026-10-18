# Rank-one factor of the ordinal latent scores' covariance

Writes `cov` as \\l l^T + I\\ when `cov - I` has rank at most one: a
single ordinal random effect, as seen through the patient's ordinal
rows. Returns `NULL` (use `pmvnorm()`) otherwise, including when
`cov - I` is not positive semi-definite (the latent residual variance is
then not 1).

## Usage

``` r
latent_factor(cov)
```

## Value

The vector \\l\\ (all zero for rank zero), or `NULL`.
