# Conditional distribution of Y\|T for a mixed continuous/ordinal (Gaussian-copula) joint model, if no competing risk – biomarker-value prediction

Copula-aware counterpart to
[`conditionalYTBio()`](https://wenhaoli18.github.io/BJM/reference/conditionalYTBio.md),
used by
[`dynamicPredictionBio()`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md)
whenever `long_fit_all$biomarker_type` contains at least one `"ordinal"`
biomarker (see
[`longitudinalSubCopula()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSubCopula.md)).
Both a continuous and an **ordinal** `bio_i` target are supported: for
an ordinal `bio_i`, `Y_all` is expected to already be the vector of
candidate *category labels* (the fitted factor's
[`levels()`](https://rdrr.io/r/base/levels.html), in threshold order)
rather than a numeric grid – see
[`compute_bio_marker_step()`](https://wenhaoli18.github.io/BJM/reference/compute_bio_marker_step.md),
which builds that vector and translates the result back into integer
category codes for its caller. No branching is actually needed here:
[`build_conditional_design_copula()`](https://wenhaoli18.github.io/BJM/reference/build_conditional_design_copula.md)
already looks up each row of biomarker `bio_i` (historical **and** the
candidate row assigned below) generically via
`long_fit_all$biomarker_type[bio_i]`, bracketing an ordinal row's latent
score between its category's cumulative-link thresholds instead of
matching it to an exact value – exactly the same mechanism already used
for every *other* ordinal biomarker's observed history in the joint
density. Assigning a candidate label into `data_it_Y`'s ordinal factor
column (below) therefore evaluates the same
[`mixed_density_prob_copula()`](https://wenhaoli18.github.io/BJM/reference/mixed_density_prob_copula.md)
box probability as any other ordinal row would, with no separate code
path required.

Unlike
[`conditionalYTBio()`](https://wenhaoli18.github.io/BJM/reference/conditionalYTBio.md)
– which can evaluate
[`mvtnorm::dmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/Mvnorm.html) once
per `l_i` across every `Y_all` candidate in one vectorized call, because
the multivariate-normal density does not need re-normalizing per
candidate – here each candidate value in `Y_all` changes the conditional
block's mean/covariance (see
[`mixed_density_prob_copula()`](https://wenhaoli18.github.io/BJM/reference/mixed_density_prob_copula.md))
and therefore requires its own
[`mvtnorm::pmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/pmvnorm.html)
evaluation. This makes the copula path in this function
`O(patients * length(l_i) * length(Y_all))` Monte-Carlo box probability
evaluations, rather than one `dmvnorm()` call per patient/`l_i` – keep
`bandcount3` modest for mixed fits with a continuous `bio_i` (an ordinal
`bio_i`'s `Y_all` length is fixed at its category count, not controlled
by `bandcount3` at all).

## Usage

``` r
conditionalYTBioCopula(
  Y_all,
  time_new,
  bio_i,
  data_predict_all,
  long_fit_all,
  l_i,
  survival_variable,
  time_variable,
  survival_variable_all,
  survival_trans_function
)
```

## Arguments

- time_new:

  Prediction time add horizon

- bio_i:

  Biomarker used to do prediction

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
