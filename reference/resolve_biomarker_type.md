# Resolve the per-biomarker type used by `longitudinalSub()`

If `biomarker_type` is supplied explicitly, it always wins (after
validation). Otherwise each biomarker's type is auto-detected from its
own response column in `data_fit_all`: a `factor` (or `ordered factor`)
response is treated as `"ordinal"`; anything else is treated as
`"continuous"`. This mirrors exactly the priority rule requested for
[`longitudinalSub`](https://liwh0904.github.io/BJM/reference/longitudinalSub.md):
explicit `biomarker_type` overrides auto-detection, and auto-detection
is only consulted when `biomarker_type` is `NULL`.

## Usage

``` r
resolve_biomarker_type(biomarker_type, data_fit_all, long_sub_fixed, M)
```
