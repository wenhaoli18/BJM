# Can a model frame's survival-time columns be overwritten in place?

[`conditionalYDT()`](https://wenhaoli18.github.io/BJM/reference/conditionalYDT.md)
and
[`conditionalYDTBio()`](https://wenhaoli18.github.io/BJM/reference/conditionalYDTBio.md)
build each biomarker's model frame once per patient and, for every
integration grid point, overwrite its survival-time columns in place
instead of rebuilding it. That is only valid when the survival time (and
every `survival_variable_all` column) enters the model frame as a bare
column: a model-frame column such as `log(years)`, `I(years^2)` or
`poly(years, 2)` is named by its expression, so it was never found and
kept its value at the first grid point for the whole integral. This
reports whether every model-frame variable that refers to one of
`time_vars` is one of those bare columns.

## Usage

``` r
time_columns_bare(terms_model, time_vars)
```

## Arguments

- terms_model:

  The fitted model's `terms`.

- time_vars:

  Names of the survival-time column and its user-supplied
  transformations.

## Value

`TRUE` if in-place updating is valid.
