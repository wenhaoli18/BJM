# Compare two prediction results' plain numeric-vector fields

Shared comparison logic for
[`checkBandcountConvergence()`](https://liwh0904.github.io/BJM/reference/checkBandcountConvergence.md)
and the `"auto"` bandcount support in
[`predictRisk()`](https://liwh0904.github.io/BJM/reference/predictRisk.md)/
[`dynamicPredictionBio()`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBio.md):
only plain numeric vectors (no `dim`) that have the same length in both
results are compared. This naturally skips fields whose *size* is itself
controlled by the bandcount being varied (e.g.
[`dynamicPredictionBio()`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBio.md)'s
`Y_density` matrix and `Y_all` grid, whose resolution is exactly what
`bandcount3` sets), while still comparing the actual per-patient
estimates derived from them (`risk_prob_1`/`risk_prob_2`, `Y_predict`).

## Usage

``` r
max_relative_diff(result_a, result_b)
```

## Value

A list with `max` (the largest relative change across all comparable
fields, or `NA` if none were comparable) and `by_field` (a named numeric
vector, one entry per comparable field).
