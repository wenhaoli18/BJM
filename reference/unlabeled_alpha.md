# Marginal-likelihood estimate of the log baseline-hazard shift

Computes each patient's log density \\\log \sum_D f(Y \mid T, D) P(D
\mid T)\\ on a grid of event times once – it does not depend on
\\\alpha\\ – and then maximizes the marginal log likelihood over
\\\alpha\\, which only moves the event-time probabilities on the grid.

## Usage

``` r
unlabeled_alpha(
  source_fit,
  data_predict_all,
  long_fit_all,
  time_variable,
  survival_variable_all,
  survival_trans_function,
  landmark,
  n_grid,
  upper,
  alpha_range,
  cap_event_time
)
```
