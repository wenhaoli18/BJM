# Outcome variables of the survival sub-model

The survival-time variable of `survival_fit_all` and, with competing
risks, its event-type variable (the response of `form_conditional_cr`).
Both are unknown for a patient still at risk at `prediction_time`; the
prediction functions integrate over them.

## Usage

``` r
outcome_variables(survival_fit_all)
```

## Arguments

- survival_fit_all:

  Output of
  [`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md).

## Value

A character vector of variable names.
