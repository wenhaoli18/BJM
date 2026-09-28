# Normalize a risk-probability ratio into a valid probability

Shared helper for `predictRisk` and `dynamicPredictionBio`: divides
summed predicted-event mass by summed total mass and clamps the result
to `[0, 1]`.

## Usage

``` r
clamp_risk_prob(numerator_sum, denominator_sum)
```

## Value

A numeric vector of risk probabilities in `[0, 1]`.
