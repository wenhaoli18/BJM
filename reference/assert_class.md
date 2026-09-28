# Assert that an object is the output of a specific BJM fitting function

Shared input-validation helper: checks the S3 class tag attached by
[`survivalSub()`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md)/[`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md),
so passing the wrong object (or the arguments in the wrong order) fails
immediately with a clear message instead of deep inside the prediction
code.

## Usage

``` r
assert_class(x, expected_class, arg_name, expected_source)
```
