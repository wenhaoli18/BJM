# Conditional distribution of Y\|D,T for a mixed continuous/ordinal (Gaussian-copula) joint model, if with competing risk

Copula-aware counterpart to
[`conditionalYDT()`](https://wenhaoli18.github.io/BJM/reference/conditionalYDT.md),
used by
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
whenever `long_fit_all$biomarker_type` contains at least one `"ordinal"`
biomarker (see
[`longitudinalSubCopula()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSubCopula.md)).
See
[`conditionalYTCopula()`](https://wenhaoli18.github.io/BJM/reference/conditionalYTCopula.md)
for the mixed continuous/ordinal density/probability construction shared
by both event-type branches (`w0`/`w1`) here.

## Usage

``` r
conditionalYDTCopula(
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
