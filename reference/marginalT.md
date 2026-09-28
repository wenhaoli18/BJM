# Marginal distribution of T

Marginal distribution of T

## Usage

``` r
marginalT(data_predict_all, long_fit_all, survival_fit_all, l_i, upper_bound)
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

- upper_bound:

  Upper limit for integration. To manage the hazard function,
  extrapolation is required. The `upper_bound` parameter specifies the
  upper time points at which extrapolation is performed.

## Value

Marginal density probability of the survival variable T is represented
as a probability matrix. In this matrix, the rows (l_i) are aligned with
specific time points, while the columns correspond to individual
patients. Each entry in the matrix denotes the marginal density
probability of survival for a given patient at a particular time point.
