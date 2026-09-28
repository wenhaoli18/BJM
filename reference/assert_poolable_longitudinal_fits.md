# Assert that a list of longitudinalSub.BJM fits is poolable with Rubin's rules

Shared input-validation helper for
[`poolLongitudinalSub`](https://wenhaoli18.github.io/BJM/reference/poolLongitudinalSub.md):
checks that `long_fit_all_list` is a list of at least two
`longitudinalSub.BJM` objects (the S3 class
[`longitudinalSub`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
attaches to its return value), and that every one of them was fit for
the same set of biomarkers with the same fixed-effect coefficients –
otherwise "the mean of the estimates" would be averaging unrelated
quantities across completions.

## Usage

``` r
assert_poolable_longitudinal_fits(
  long_fit_all_list,
  arg_name = "long_fit_all_list"
)
```
