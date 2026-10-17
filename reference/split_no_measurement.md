# Separate the patients `prepare_risk_data()` found without measurements

A patient none of whose rows, in any biomarker's data, has a complete
measurement
([`prepare_risk_data()`](https://wenhaoli18.github.io/BJM/reference/prepare_risk_data.md)
has already dropped the incomplete rows of every other patient) is
predicted from the survival sub-model alone.

## Usage

``` r
split_no_measurement(
  data_predict_all,
  long_fit_all,
  survival_variable,
  survival_variable_all
)
```

## Arguments

- survival_variable:

  The outcome variables of the survival sub-model
  ([`outcome_variables()`](https://wenhaoli18.github.io/BJM/reference/outcome_variables.md)),
  which are not checked, as in
  [`drop_missing_longitudinal()`](https://wenhaoli18.github.io/BJM/reference/drop_missing_longitudinal.md).

## Value

A list: `measured`, `data_predict_all` without these patients, and
`no_measurement`, a list with one row per such patient in each element
(as
[`marginalT()`](https://wenhaoli18.github.io/BJM/reference/marginalT.md)
and friends read covariates from), or `NULL` if there is none.
