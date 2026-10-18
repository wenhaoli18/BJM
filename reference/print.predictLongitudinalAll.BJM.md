# Print method for `predictLongitudinalAll.BJM` objects

Called when you print the result of
[`predictLongitudinal()`](https://wenhaoli18.github.io/BJM/reference/predictLongitudinal.md)
with several biomarkers (`bio_i` naming more than one, or `NULL`):
prints each biomarker's prediction in turn, under its name.

## Usage

``` r
# S3 method for class 'predictLongitudinalAll.BJM'
print(x, ...)
```

## Arguments

- x:

  A `predictLongitudinalAll.BJM` object.

- ...:

  Passed to
  [`print.predictLongitudinal.BJM()`](https://wenhaoli18.github.io/BJM/reference/print.predictLongitudinal.BJM.md).

## Value

Invisibly returns `x`.
