# Auto-select "auto" bandcount3 for a single biomarker's per-marker step

[`dynamicPredictionBioAll()`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBioAll.md)'s
counterpart to
[`auto_tune_bandcount()`](https://liwh0904.github.io/BJM/reference/auto_tune_bandcount.md):
`bandcount3` only controls
[`compute_bio_marker_step()`](https://liwh0904.github.io/BJM/reference/compute_bio_marker_step.md)'s
own candidate-value grid (`Y_all`), not the shared step, so it is tuned
per biomarker by re-calling
[`compute_bio_marker_step()`](https://liwh0904.github.io/BJM/reference/compute_bio_marker_step.md)
directly against an already-computed `shared` object (from
[`compute_bio_shared_step()`](https://liwh0904.github.io/BJM/reference/compute_bio_shared_step.md)),
rather than re-running the whole
[`dynamicPredictionBio()`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBio.md)
pipeline – which would recompute the shared denominator on every
doubling round, for every biomarker, exactly the redundant work
[`dynamicPredictionBioAll()`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBioAll.md)
is meant to avoid. Same doubling-until-stable check as
[`auto_tune_bandcount()`](https://liwh0904.github.io/BJM/reference/auto_tune_bandcount.md)
(compares `Y_predict` via
[`max_relative_diff()`](https://liwh0904.github.io/BJM/reference/max_relative_diff.md),
capped at `max_rounds` doublings, warns instead of erroring if not
converged by then).

## Usage

``` r
auto_tune_marker_bandcount3(
  shared,
  bio_i,
  long_fit_all,
  survival_fit_all,
  prediction_time,
  horizon,
  time_variable,
  survival_variable_all,
  survival_trans_function,
  tol = 0.01,
  max_rounds = 2,
  multiplier = 2
)
```

## Arguments

- shared:

  Output of
  [`compute_bio_shared_step()`](https://liwh0904.github.io/BJM/reference/compute_bio_shared_step.md).

- bio_i:

  See
  [`dynamicPredictionBio`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBio.md).

- long_fit_all:

  See
  [`dynamicPredictionBio`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBio.md).

- survival_fit_all:

  See
  [`dynamicPredictionBio`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBio.md).

- prediction_time:

  See
  [`dynamicPredictionBio`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBio.md).

- horizon:

  See
  [`dynamicPredictionBio`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBio.md).

- time_variable:

  See
  [`dynamicPredictionBio`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBio.md).

- survival_variable_all:

  See
  [`dynamicPredictionBio`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBio.md).

- survival_trans_function:

  See
  [`dynamicPredictionBio`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBio.md).

## Value

A list with `result`
([`compute_bio_marker_step()`](https://liwh0904.github.io/BJM/reference/compute_bio_marker_step.md)'s
return value at the resolved `bandcount3`) and `bandcount3` (the
resolved numeric value).
