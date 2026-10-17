# Default covariate update for the episodes after the open one

Adds `k` to the episode counter, but never past `cap`: the largest
episode number in the data the survival sub-model was fit on. Beyond it
the model's effect of the counter is an extrapolation, and a hazard that
rises with each earlier recurrence, run out to episode 30, multiplies
the hazard many hundredfold, so paths drag on through ever shorter
episodes and the counts are inflated.

## Usage

``` r
default_next_covariates(episode_variable, cap = Inf)
```
