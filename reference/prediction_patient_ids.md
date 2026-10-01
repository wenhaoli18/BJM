# Patient ids, in the order predictions are returned

The density helpers
([`marginalT()`](https://wenhaoli18.github.io/BJM/reference/marginalT.md),
[`conditionalYT()`](https://wenhaoli18.github.io/BJM/reference/conditionalYT.md),
...) all produce one column per patient, in order of first appearance in
`data_predict_all[[1]]`; this returns those ids, as character, to name
the returned predictions with.

## Usage

``` r
prediction_patient_ids(data_predict_all, long_fit_all)
```

## Value

A character vector.
