# Bio_i-specific computation for dynamicPredictionBio()

Internal helper for the part of
[`dynamicPredictionBio()`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBio.md)'s
pipeline that does depend on which biomarker is being predicted:
building the candidate-value grid `Y_all`, evaluating the numerator
conditional density on it, and assembling/normalizing `Y_density` and
`Y_predict`. Takes the output of
[`compute_bio_shared_step()`](https://liwh0904.github.io/BJM/reference/compute_bio_shared_step.md)
so the bio_i-independent pieces are not recomputed.

## Usage

``` r
compute_bio_marker_step(
  shared,
  bio_i,
  long_fit_all,
  survival_fit_all,
  prediction_time,
  horizon,
  time_variable,
  survival_variable_all,
  survival_trans_function,
  bandcount3
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

- bandcount3:

  See
  [`dynamicPredictionBio`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBio.md).

## Value

A list with `Y_predict`, `Y_density`, `Y_all` (matching
[`dynamicPredictionBio()`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBio.md)'s
return value fields). For an **ordinal** `bio_i`, `Y_all`/`Y_predict`
hold integer category codes (`1:K`, in threshold order) rather than a
numeric grid – see Details below – and `Y_all` additionally carries a
`"category_labels"` attribute with the matching level-label strings.
