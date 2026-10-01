# Set a patient's survival-time columns to one integration grid point

Sets the survival-time column and each `survival_variable_all` column
(through its `survival_trans_function`) to their values at event time
`l`. With `only_existing = TRUE` (for a model frame updated in place,
see
[`time_columns_bare()`](https://wenhaoli18.github.io/BJM/reference/time_columns_bare.md))
columns not already present are left out.

## Usage

``` r
set_survival_columns(
  df,
  survival_variable,
  l,
  survival_variable_all,
  survival_trans_function,
  only_existing = FALSE
)
```

## Value

`df`, updated.
