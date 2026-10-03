# Draw from a truncated multivariate normal distribution

`n` draws of \\N(\mu, V)\\ restricted to the box `(lower, upper]`: by
rejection from batches of untruncated proposals, and, for any draws
still missing once `max_proposals` proposals have been spent (a box of
small probability), by Gibbs sampling, one coordinate at a time from its
univariate truncated normal conditional distribution, after `burn_in`
sweeps and keeping every `thin`-th one.

## Usage

``` r
rtmvnorm_box(
  n,
  mean,
  sigma,
  lower,
  upper,
  max_proposals = 20000,
  burn_in = 50,
  thin = 5
)
```

## Value

A matrix with one column per draw.
