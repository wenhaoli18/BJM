# A biomarker's model frame at one integration grid point

Overwrites the survival-time columns of the per-patient template `mf` in
place when `in_place` (see
[`time_columns_bare()`](https://wenhaoli18.github.io/BJM/reference/time_columns_bare.md));
otherwise rebuilds the model frame from `data` with those columns set,
so that transformed terms such as `log(years)` are re-evaluated at `l`.

## Usage

``` r
survival_model_frame_at(
  mf,
  data,
  in_place,
  terms_model,
  xlev,
  survival_variable,
  l,
  survival_variable_all,
  survival_trans_function
)
```

## Value

The model frame.
