# Bio_i-independent shared computation for dynamicPredictionBio()

Internal helper factoring out the part of
[`dynamicPredictionBio()`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBio.md)'s
pipeline that does not depend on which biomarker (`bio_i`) is being
predicted: restricting to at-risk patients, building the survival-side
integration grid out to `upper_bound`, dispatching to the
Gaussian-copula conditional-density variants when needed (see
[`predictRisk()`](https://liwh0904.github.io/BJM/reference/predictRisk.md)),
and evaluating the denominator conditional density
(`conditionalYT`/`conditionalYDT`).
[`dynamicPredictionBioAll()`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBioAll.md)
computes this once and reuses it across every requested biomarker,
instead of recomputing it once per biomarker (including, under the
copula path, its
[`mvtnorm::pmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/pmvnorm.html)
calls) the way calling
[`dynamicPredictionBio()`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBio.md)
once per biomarker would.

## Usage

``` r
compute_bio_shared_step(
  data_predict_all,
  long_fit_all,
  survival_fit_all,
  prediction_time,
  time_variable,
  survival_variable_all,
  survival_trans_function,
  bandcount2
)
```

## Arguments

- data_predict_all:

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

- time_variable:

  See
  [`dynamicPredictionBio`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBio.md).

- survival_variable_all:

  See
  [`dynamicPredictionBio`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBio.md).

- survival_trans_function:

  See
  [`dynamicPredictionBio`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBio.md).

- bandcount2:

  See
  [`dynamicPredictionBio`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBio.md).

## Value

A list with the pieces
[`compute_bio_marker_step()`](https://liwh0904.github.io/BJM/reference/compute_bio_marker_step.md)
needs: `data_predict_all` (at-risk-filtered), `survival_variable`,
`predict.time.infinity`, `S_T_all_infinity`, `has_cr`,
`D_T_all_infinity` (`NULL` if no competing risk), `f_y_D_all_infinity`,
and the dispatched `conditionalYTBio_fun`/`conditionalYDTBio_fun`
functions.
