# Standard deviation of each patient's predictive distribution

Helper for
[`max_relative_diff()`](https://wenhaoli18.github.io/BJM/reference/max_relative_diff.md):
from a
[`dynamicPredictionBio()`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md)-style
result (`Y_density`, one column per patient, tabulated on `Y_all`), the
standard deviation of each patient's predicted biomarker distribution.
`0` if the result has no density.

## Usage

``` r
predictive_sd(result)
```

## Arguments

- result:

  A prediction result.

## Value

A numeric vector, one per patient, or `0`.
