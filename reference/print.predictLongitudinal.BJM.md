# Print method for `predictLongitudinal.BJM` objects

Automatically called when you type the result of
[`predictLongitudinal()`](https://wenhaoli18.github.io/BJM/reference/predictLongitudinal.md)
(with a single `bio_i`) at the console.

## Usage

``` r
# S3 method for class 'predictLongitudinal.BJM'
print(
  x,
  bio_i = NULL,
  long_fit_all = NULL,
  prediction_time = NULL,
  horizon = NULL,
  subject_ids = NULL,
  digits = 4,
  ...
)
```

## Arguments

- x:

  A `predictLongitudinal.BJM` object.

- bio_i:

  Biomarker index (for label lookup). Default `NULL`.

- long_fit_all:

  `longitudinalSub.BJM` object for name lookup.

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
