# Filter longitudinal data down to a single patient, substituting the candidate biomarker prediction value

Shared helper for `conditionalYTBio` and `conditionalYDTBio`: extracts
each biomarker's rows for one patient ID, and for the biomarker being
predicted (`bio_i`), appends a row at `time_new` with each candidate
value in `Y_all` substituted in turn.

## Usage

``` r
select_patient_longitudinal_data_bio(
  data.long,
  num,
  num_i,
  n_longitudinal,
  time_variable,
  bio_i,
  time_new,
  Y_all,
  long_fit_all
)
```

## Value

A list with `rep_num_i_list`, `data_num_i_list`, and `Y_select_all` (a
matrix of the substituted biomarker values, one column per element of
`Y_all`), or `NULL` if any biomarker has zero rows for this patient
(caller should skip the patient in that case).
