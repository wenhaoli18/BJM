# Print method for `predictRisk.BJM` objects

Automatically called when you type the result of
[`predictRisk()`](https://liwh0904.github.io/BJM/reference/predictRisk.md)
at the console.

## Usage

``` r
# S3 method for class 'predictRisk.BJM'
print(
  x,
  prediction_time = NULL,
  horizon = NULL,
  subject_ids = NULL,
  digits = 4,
  ...
)
```

## Arguments

- x:

  A `predictRisk.BJM` object.

- prediction_time:

  Landmark time (for display). Default `NULL`.

- horizon:

  Prediction horizon (for display). Default `NULL`.

- subject_ids:

  Optional subject ID labels.

- digits:

  Decimal places. Default 4.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `x`.
