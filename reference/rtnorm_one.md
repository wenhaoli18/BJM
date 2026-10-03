# Draw from a univariate truncated normal distribution

By inversion, on the side of the mean where the truncation interval lies
(so a box far in the upper tail does not lose precision to
[`pnorm()`](https://rdrr.io/r/stats/Normal.html) rounding to 1).

## Usage

``` r
rtnorm_one(mean, sd, lower, upper)
```
