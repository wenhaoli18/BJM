# Factor levels of the data a biomarker's sub-model was fit on

Used by
[`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
to cache, per biomarker, the factor levels seen by `lme()` – computed
from the full fitting data before it is restricted to the subjects
retained in the joint fit, which may lack some levels.

## Usage

``` r
training_xlevels(terms_model, data, fixed_formula)
```

## Arguments

- terms_model:

  The fitted model's `terms`.

- data:

  The data the model was fit on.

- fixed_formula:

  The fixed-effects formula (for its variables).

## Value

A named list of factor levels, as from
[`.getXlevels()`](https://rdrr.io/r/stats/checkMFClasses.html).
