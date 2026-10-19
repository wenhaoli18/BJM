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
  [`poolLongitudinalSub`](https://wenhaoli18.github.io/BJM/reference/poolLongitudinalSub.md).

- digits:

  Number of significant digits. Default is 4.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `x`.

## Examples

``` r
# \donttest{
if (requireNamespace("torch", quietly = TRUE) && torch::torch_is_installed()) {
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
#> 
#> Call:
#> poolLongitudinalSub(m = 5 completed-data fits, M = 2 longitudinal outcomes)
#> 
#> =================================================================
#>  Rubin's-Rules-Pooled Longitudinal Sub-model(s)
#> -----------------------------------------------------------------
#>  [1] long1
#> 
#>                Value       SE    DF       t  p-value    
#> (Intercept)  2.63947  0.37760 697.3   6.990 6.43e-12 ***
#> year         0.26795  0.01842 611.0  14.549  < 2e-16 ***
#> age         -0.01296  0.00572 162.6  -2.265   0.0248 *  
#> sex          0.04997  0.16523 162.8   0.302   0.7627    
#> years       -0.22939  0.02095 161.5 -10.949  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>   Fraction of missing information (FMI), by coefficient:
#>   (Intercept) = 0.006,  year = 0.029,  age = 0.014,  sex = 0.013,  years = 0.019 
#> 
#> -----------------------------------------------------------------
#>  [2] long2
#> 
#>                 Value        SE    DF       t p-value    
#> (Intercept)  3.531148  0.151814 598.3  23.260 < 2e-16 ***
#> year        -0.141220  0.010241 599.2 -13.789 < 2e-16 ***
#> age         -0.006819  0.002334 146.5  -2.922 0.00403 ** 
#> sex         -0.163901  0.066921 150.6  -2.449 0.01546 *  
#> years        0.085807  0.008517 151.4  10.075 < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>   Fraction of missing information (FMI), by coefficient:
#>   (Intercept) = 0.032,  year = 0.031,  age = 0.054,  sex = 0.047,  years = 0.045 
#> 
#> =================================================================
# }
```
