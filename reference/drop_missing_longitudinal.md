# Drop prediction rows with a missing biomarker value or covariate

Shared helper for `predictRisk`, `dynamicPredictionBio`, and
`dynamicPredictionBioAll`: a patient's conditional density only involves
the measurements actually observed, so for each biomarker, rows of its
`data_predict_all` element with a missing response, or a missing
covariate/time/ID used by that biomarker's
`long_sub_fixed`/`long_sub_random` formulas, are removed. Previously
such a row put an `NA` into the stacked outcome vector (or misaligned
the design matrix) and made that patient's whole prediction `NA`. The
outcome variables (`survival_variable`: the survival time and, with
competing risks, the event type – see
[`outcome_variables()`](https://wenhaoli18.github.io/BJM/reference/outcome_variables.md))
and `survival_variable_all` are not checked: they are overwritten with
each integration grid point / event type before use, and are unknown
(typically `NA`) for a patient still at risk. The event type used to be
checked, so with the event type in `long_sub_fixed` every at-risk
patient was dropped. A patient left with no rows for some biomarker
cannot be predicted; a warning names them, and it is an error if that
leaves no patient at all.

## Usage

``` r
drop_missing_longitudinal(
  data_predict_all,
  long_fit_all,
  survival_variable,
  survival_variable_all
)
```

## Value

`data_predict_all`, filtered per element.
