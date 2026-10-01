# Random-effects design matrix for one patient, across all biomarkers

Shared helper for
[`build_conditional_design()`](https://wenhaoli18.github.io/BJM/reference/build_conditional_design.md),
[`build_conditional_design_copula()`](https://wenhaoli18.github.io/BJM/reference/build_conditional_design_copula.md),
and
[`process_variance()`](https://wenhaoli18.github.io/BJM/reference/process_variance.md):
builds each biomarker's random-effects design from the left-hand side of
its own `long_sub_random` formula (e.g. `~ 1`, `~ year`,
`~ year + I(year^2)`) and stacks them block-diagonally, matching the
block-diagonal `Sigma_fit` built by
[`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md).
Previously this assumed every biomarker had either a random intercept
only, or a random intercept plus a slope on `time_variable`, so a mix of
the two (or any other random-effects formula) failed with
"non-conformable arguments".

## Usage

``` r
random_effects_design(data_num_i_list, long_sub_random, Sigma)
```

## Arguments

- data_num_i_list:

  One patient's data, one data frame per biomarker.

- long_sub_random:

  The list of random-effects formulas.

- Sigma:

  The random-effects covariance matrix (`Sigma_fit`).

## Value

A matrix with one row per observation (stacked across biomarkers) and
one column per random effect.
