# Conditional distribution of Y\|D,T for a mixed continuous/ordinal (Gaussian-copula) joint model, if with competing risk – biomarker-value prediction

Copula-aware counterpart to
[`conditionalYDTBio()`](https://wenhaoli18.github.io/BJM/reference/conditionalYDTBio.md),
used by
[`dynamicPredictionBio()`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md)
whenever `long_fit_all$biomarker_type` contains at least one `"ordinal"`
biomarker and a competing-risk `survival_fit_all` is supplied. Both a
continuous and an **ordinal** `bio_i` target are supported: for an
ordinal `bio_i`, `Y_all` is expected to already be the vector of
candidate *category labels* (the fitted factor's
[`levels()`](https://rdrr.io/r/base/levels.html), in threshold order)
rather than a numeric grid – see
[`compute_bio_marker_step()`](https://wenhaoli18.github.io/BJM/reference/compute_bio_marker_step.md),
which builds that vector and translates the result back into integer
category codes for its caller. See
[`conditionalYTBioCopula()`](https://wenhaoli18.github.io/BJM/reference/conditionalYTBioCopula.md)
for why no branching is actually needed here –
[`build_conditional_design_copula()`](https://wenhaoli18.github.io/BJM/reference/build_conditional_design_copula.md)
already looks up each row of biomarker `bio_i` generically via
`long_fit_all$biomarker_type[bio_i]`, whether that row is historical or
the candidate row assigned below – and for why this function evaluates
[`mvtnorm::pmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/pmvnorm.html)
once per patient/`l_i`/ `Y_all` candidate/event-type branch, rather than
one vectorized `dmvnorm()` call per patient/`l_i` as
[`conditionalYDTBio()`](https://wenhaoli18.github.io/BJM/reference/conditionalYDTBio.md)
does.

## Usage

``` r
conditionalYDTBioCopula(
  Y_all,
  time_new,
  bio_i,
  data_predict_all,
  long_fit_all,
  survival_fit_all,
  l_i,
  survival_variable,
  time_variable,
  survival_variable_all,
  survival_trans_function
)
```

## Arguments

- data_predict_all:

  This involves a collection of `data.frame` objects for dynamic
  prediction, each corresponding to a distinct longitudinal outcome.
  These data frames should contain the variables specified in
  `long_sub_fixed` and `long_sub_random`. Utilizing a list structure
  allows for the incorporation of multiple longitudinal outcomes, each
  potentially following different measurement protocols. In instances
  where all longitudinal outcomes are recorded at identical time points
  across patients, a singular `data.frame` object may be used in a
  `list`. It is presumed that each data frame is structured in a long
  format.

- long_fit_all:

  Outputs from the model fitting process using the `nlme` package,
  encompassing the results and parameters obtained from the analysis.

- l_i:

  A vector of time points to calculate the conditional probability.

- survival_variable:

  Time-to-event outcomes variable name.

- time_variable:

  The name of time variable in linear mixed model.

- survival_variable_all:

  The name of the transformed time-to-event outcomes variable.

- survival_trans_function:

  The transformation function used for time-to-event outcomes, in the
  order of `survival_variable_all`.
