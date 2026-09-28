# Conditional distribution of Y\|T for a mixed continuous/ordinal (Gaussian-copula) joint model, if no competing risk

Copula-aware counterpart to
[`conditionalYT()`](https://liwh0904.github.io/BJM/reference/conditionalYT.md),
used by
[`predictRisk()`](https://liwh0904.github.io/BJM/reference/predictRisk.md)
whenever `long_fit_all$biomarker_type` contains at least one `"ordinal"`
biomarker (see
[`longitudinalSubCopula()`](https://liwh0904.github.io/BJM/reference/longitudinalSubCopula.md)).
Continuous biomarkers contribute an exact Gaussian density factor,
exactly as
[`conditionalYT()`](https://liwh0904.github.io/BJM/reference/conditionalYT.md)
computes for every biomarker; ordinal biomarkers instead contribute a
Gaussian-copula box probability (via
[`mixed_density_prob_copula()`](https://liwh0904.github.io/BJM/reference/mixed_density_prob_copula.md)),
since only the cumulative-link category – not the exact underlying
latent score – is observed for them. All-continuous fits are unaffected:
they are still routed to
[`conditionalYT()`](https://liwh0904.github.io/BJM/reference/conditionalYT.md)
by
[`predictRisk()`](https://liwh0904.github.io/BJM/reference/predictRisk.md),
this function is never called for them, and
[`mixed_density_prob_copula()`](https://liwh0904.github.io/BJM/reference/mixed_density_prob_copula.md)
reduces to (a constant multiple of)
[`conditionalYT()`](https://liwh0904.github.io/BJM/reference/conditionalYT.md)'s
own computation when there are no ordinal markers – see that function's
documentation.

## Usage

``` r
conditionalYTCopula(
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
