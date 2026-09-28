# Build the longitudinal outcome matrix for all candidate biomarker values

Shared helper for `conditionalYTBio` and `conditionalYDTBio`: for each
candidate value in `Y_all`, assembles the row of observed longitudinal
outcomes across biomarkers, substituting the candidate value for the
biomarker being predicted.

## Usage

``` r
build_longitudinal_matrix_bio(
  data_num_i_list,
  lfit,
  bio_i,
  Y_select_all,
  n_longitudinal,
  Y_all
)
```

## Value

A matrix with `length(Y_all)` rows, one per candidate value.
