# Check the event-time transformations against the training data

The longitudinal sub-models were fit on data whose
`survival_variable_all` columns the user computed from the event time;
at prediction time the same columns are recomputed with
`survival_trans_function`. The two are passed as separate lists, so a
function out of order (or a different one) used to give wrong
predictions silently. This applies each function to the event times
stored with every fit (`lme()`'s `data`, `clmm()`'s model frame) and
stops if the result differs from the column the model was fit on.

## Usage

``` r
check_survival_trans_training(
  long_fit_all,
  survival_fit_all,
  survival_variable_all,
  survival_trans_function
)
```

## Arguments

- long_fit_all:

  Output of
  [`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md).

- survival_fit_all:

  Output of
  [`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md).

- survival_variable_all, survival_trans_function:

  As in
  [`predictRisk`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md).

## Value

`TRUE`, invisibly; an error on a mismatch.
