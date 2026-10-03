# A patient's history rows usable for conditioning

The rows of one biomarker's data with the biomarker and every covariate
of its sub-model observed (the outcome columns, which are set at each
candidate event time, excepted). May have no rows.

## Usage

``` r
complete_history_rows(d, i, long_fit_all, survival_variable, other_outcomes)
```
