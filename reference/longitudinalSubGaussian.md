# Fit a multivariate longitudinal sub-model (all-continuous biomarkers)

Internal workhorse for
[`longitudinalSub`](https://liwh0904.github.io/BJM/reference/longitudinalSub.md)
when every biomarker is continuous. This is the original
`longitudinalSub` implementation, kept verbatim and called with the
caller's pristine, unmodified arguments, so that the all-continuous case
behaves exactly as it always has – see
[`longitudinalSub`](https://liwh0904.github.io/BJM/reference/longitudinalSub.md)
for the argument documentation, and
[`longitudinalSubCopula()`](https://liwh0904.github.io/BJM/reference/longitudinalSubCopula.md)
for the mixed continuous/ordinal extension used when at least one
biomarker is categorical.

## Usage

``` r
longitudinalSubGaussian(data_fit_all, long_sub_fixed, long_sub_random)
```
