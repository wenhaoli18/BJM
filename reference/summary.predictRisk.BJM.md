# Summary method for `predictRisk.BJM` objects

Like `print` but also shows mean, SD, and range of predicted risks.

## Usage

``` r
# S3 method for class 'predictRisk.BJM'
summary(
  object,
  prediction_time = NULL,
  horizon = NULL,
  subject_ids = NULL,
  digits = 4,
  ...
)
```

## Arguments

- object:

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

Invisibly returns `object`.
