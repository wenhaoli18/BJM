# One patient's stacked fixed-effects means at every grid point

[`fixed_mean_grid()`](https://wenhaoli18.github.io/BJM/reference/fixed_mean_grid.md)
for each biomarker's rows of one patient (`data_num_i_list`), stacked in
biomarker order as in
[`build_conditional_design()`](https://wenhaoli18.github.io/BJM/reference/build_conditional_design.md).

## Usage

``` r
patient_mean_grid(
  data_num_i_list,
  long_fit_all,
  l_i,
  survival_variable,
  survival_variable_all,
  survival_trans_function,
  fixed_values = list()
)
```

## Value

A matrix, one row per stacked observation and one column per grid point
in `l_i`.
