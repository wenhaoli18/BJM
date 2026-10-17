# Log densities of an empty longitudinal history

Stand-ins for
[`conditionalYT()`](https://wenhaoli18.github.io/BJM/reference/conditionalYT.md)/[`conditionalYDT()`](https://wenhaoli18.github.io/BJM/reference/conditionalYDT.md)
for patients without any measurement: \\f(Y \mid T) = 1\\, i.e. a log
density of 0 at every time (and event type), so
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
gives the survival sub-model's risk.

## Usage

``` r
conditionalYT_none(data_predict_all, long_fit_all, l_i, ...)

conditionalYDT_none(data_predict_all, long_fit_all, survival_fit_all, l_i, ...)
```

## Value

A list of one (`conditionalYT_none()`) or two (`conditionalYDT_none()`)
zero matrices, one row per time in `l_i` and one column per patient.
