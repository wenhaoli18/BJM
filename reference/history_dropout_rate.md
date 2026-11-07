# Share of rows from subjects never measured for each biomarker

The swap rate of
[`history_dropout_inputs()`](https://wenhaoli18.github.io/BJM/reference/history_dropout_inputs.md):
for each biomarker, the share of rows from subjects with no observed
value of it, at least `min_rate`, or 0 when there is no such subject.

## Usage

``` r
history_dropout_rate(ids, mask, min_rate = 0.1)
```
