# Assert that biomarker_type is a valid per-biomarker type vector

Shared input-validation helper for the `biomarker_type` argument of
[`longitudinalSub`](https://liwh0904.github.io/BJM/reference/longitudinalSub.md):
when supplied explicitly (as opposed to being auto-detected from each
biomarker's response column), it must be a character vector with exactly
one entry per longitudinal outcome, and every entry must be one of
`"continuous"` or `"ordinal"`.

## Usage

``` r
assert_biomarker_type(
  biomarker_type,
  M,
  valid_types = c("continuous", "ordinal")
)
```
