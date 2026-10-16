# Covariates of the imputation model

The right-hand-side variables of `long_sub_fixed` other than the subject
id, including the visit time `time_variable` when a formula uses it. The
visit time is also the key that aligns rows across biomarkers, and it
used to be left out of the model along with the id, so values were
imputed without knowing when the measurement was taken (pulling a
trending biomarker towards its overall level).

## Usage

``` r
imputation_covariates(covariate_vars, id)
```

## Arguments

- covariate_vars:

  List of each biomarker's right-hand-side variables.

- id:

  Name of the id column.

## Value

Character vector of covariate names.
