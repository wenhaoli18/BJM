# Vectorized draws from truncated standard-variance normals

\\N(mean, 1)\\ truncated to \\(lower, upper\]\\, by inversion on the
side of the mean where the interval lies.

## Usage

``` r
rtnorm_vec(mean, lower, upper)
```

## Arguments

- mean, lower, upper:

  Vectors of equal length.

## Value

A numeric vector of draws.
