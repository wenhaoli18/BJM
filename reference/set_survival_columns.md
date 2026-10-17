# Set a patient's survival-time columns to one integration grid point

Sets the survival-time column and each `survival_variable_all` column
(through its `survival_trans_function`) to their values at event time
`l`.

## Usage

``` r
set_survival_columns(
  df,
  survival_variable,
  l,
  survival_variable_all,
  survival_trans_function
)
```

## Value

`df`, updated.
