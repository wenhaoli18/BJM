# Summary method for `predictLongitudinal.BJM` objects

Like `print` but also shows distribution-level summaries across
subjects.

## Usage

``` r
# S3 method for class 'predictLongitudinal.BJM'
summary(
  object,
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

- object:

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

Invisibly returns `object`.
