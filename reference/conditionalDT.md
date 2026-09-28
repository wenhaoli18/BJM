# conditional distribution of D\|T

This function computes the conditional probability density function of
competing risk event type D, given the survival time T.

## Usage

``` r
conditionalDT(data_predict_all, long_fit_all, survival_fit_all, l_i)
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

- survival_fit_all:

  Results and parameters generated from the model fitting procedure,
  utilizing the `coxph` function. These outputs include the
  comprehensive findings and variables derived from the analysis.

- l_i:

  A vector of time points to calculate the conditional probability.

## Value

Probability matrices of competing risk event type D conditional on
survival outcome T.
