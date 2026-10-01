# Variables a biomarker's sub-model uses

The rows
[`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
keeps for the joint EM must be the rows `lme()` was fit on, i.e.
complete in every variable of both the fixed- and the random-effects
formula (including the id). Only the fixed-effects variables used to be
checked, so a missing value in a variable used only in `long_sub_random`
left a row in the fixed-effects design but not in the random-effects
one, and failed with "arguments imply differing number of rows".

## Usage

``` r
model_vars(fixed_formula, random_formula)
```

## Arguments

- fixed_formula, random_formula:

  The biomarker's `long_sub_fixed` and `long_sub_random` formulas.

## Value

A character vector of variable names.
