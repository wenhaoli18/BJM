# Cox model with the generation coefficient fixed

For `survivalSubExternal(tail_hr = )`: the other coefficients are
estimated with the generation term as an offset, then the model is refit
at all coefficients with no iterations, so that the result is an
ordinary `coxph` fit (its prediction and Breslow baseline then need no
special handling of the offset).

## Usage

``` r
cox_fixed_generation(free_fit, stacked, form, generation, tail_hr)
```
