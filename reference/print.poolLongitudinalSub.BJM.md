# Print method for `poolLongitudinalSub.BJM` objects

Automatically called when you type `pooled` or `print(pooled)` at the
console. Displays a JMbayes2-style formatted summary of the
Rubin's-rules-pooled longitudinal sub-model(s).

## Usage

``` r
# S3 method for class 'poolLongitudinalSub.BJM'
print(x, digits = 4, ...)
```

## Arguments

- x:

  A `poolLongitudinalSub.BJM` object returned by
  [`poolLongitudinalSub`](https://liwh0904.github.io/BJM/reference/poolLongitudinalSub.md).

- digits:

  Number of significant digits. Default is 4.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `x`.

## Examples

``` r
# \donttest{
if (requireNamespace("torch", quietly = TRUE)) {
  data(pbc3)
  data_fit_all <- pbc3[pbc3$status3 == 1, ]

  set.seed(1)
  n <- nrow(data_fit_all)
  data_fit_all$serBilir[sample.int(n, floor(0.1 * n))] <- NA
  data_fit_all$albumin[sample.int(n, floor(0.1 * n))] <- NA

  long_sub_fixed <- list(
    "long1" = serBilir ~ year + age + sex + years,
    "long2" = albumin ~ year + age + sex + years)
  long_sub_random <- list("long1" = ~ year | id, "long2" = ~ year | id)

  imputed <- imputeLongitudinal(data_fit_all, long_sub_fixed, long_sub_random,
                                 time_variable = "year", n_imputations = 5,
                                 impute = "multiple", epochs = 50, seed = 1)

  long_fit_all_list <- lapply(imputed$data_fit_all_list, function(d) {
    longitudinalSub(d, long_sub_fixed, long_sub_random)
  })
  pooled <- poolLongitudinalSub(long_fit_all_list)
  pooled   # triggers print.poolLongitudinalSub.BJM automatically
}
#> Error: Lantern is not loaded. Please use `install_torch()` to install additional dependencies.
# }
```
