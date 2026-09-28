# Filter longitudinal data down to a single patient

Shared helper for `conditionalYT` and `conditionalYDT`: extracts each
biomarker's rows for one patient ID, replicating the shared data frame
across biomarkers when only one was supplied.

## Usage

``` r
select_patient_longitudinal_data(
  data.long,
  num,
  num_i,
  n_longitudinal,
  time_variable
)
```

## Value

A list with `rep_num_i_list` and `data_num_i_list`, or `NULL` if any
biomarker has zero rows for this patient (caller should skip the patient
in that case).
