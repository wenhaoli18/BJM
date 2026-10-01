# conditional distribution of Y\|T, if no competing risk;

This function computes the conditional probability density function of
longitudinal variable Y, given the survival time T without competing
risk D.

## Usage

``` r
conditionalYT(
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

## Value

The output is a list containing probability matrices of **log**
densities (so that they neither overflow nor underflow; see
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
for how they are exponentiated). In the presence of competing risks,
this list includes two elements; otherwise, it contains only one
element. Each element within the list is a probability matrix, with the
number of rows (l_i) corresponding to specific time points and columns
representing different patients. Every matrix element represents the
conditional probability derived from the conditional distribution of
longitudinal variable Y given the survival time T without competing risk
D for a particular patient at a specific time point.
