# Pool longitudinal sub-model fits across multiple imputations with Rubin's rules

When
[`imputeLongitudinal`](https://wenhaoli18.github.io/BJM/reference/imputeLongitudinal.md)
is run with `impute = "multiple"`, it returns `n_imputations`
independently completed versions of `data_fit_all` in
`data_fit_all_list`, instead of a single completed dataset filled with
the across-draw mean (`impute = "single"`). Fitting
[`longitudinalSub`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
once per completion and combining the resulting fixed-effect estimates
with Rubin's rules (Rubin, 1987) – rather than averaging into one
completed dataset up front – is the standard way multiple imputation
propagates the extra uncertainty from not knowing the true missing
values into the final standard errors; `impute = "single"` does not do
this (see the "Missing data imputation" article on the package website).
`poolLongitudinalSub()` performs that combination step.

## Usage

``` r
poolLongitudinalSub(long_fit_all_list)
```

## Arguments

- long_fit_all_list:

  A list of two or more `longitudinalSub.BJM` objects: the result of
  calling
  [`longitudinalSub`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)
  once on each completed dataset in
  `imputeLongitudinal(..., impute = "multiple")$data_fit_all_list`,
  using the same `long_sub_fixed`/ `long_sub_random` for every
  completion.

## Value

An object of class `"poolLongitudinalSub.BJM"`, a named list with
elements:

- pooled:

  A list, one `data.frame` per longitudinal outcome, with one row per
  fixed-effect coefficient and columns `estimate`, `std_error`, `df`,
  `statistic`, `p_value`, `riv` (relative increase in variance due to
  missingness), and `fmi` (fraction of missing information).

- m:

  The number of completed-data fits pooled.

- long_sub_fixed:

  The `long_sub_fixed` used by every fit (taken from the first element
  of `long_fit_all_list`).

- long_fit_all:

  A `longitudinalSub.BJM` object usable for prediction (e.g. as
  `long_fit_all` in
  [`predictRisk`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)):
  the first fit, with each biomarker's fixed effects replaced by their
  pooled estimates, its residual variance by the average across
  completions, and `Sigma_fit` by the average of the `m` fits'
  `Sigma_fit` (Rubin's rules point estimates). Averaging the `m` fits'
  predictions instead is a common alternative.

## Details

For each fixed-effect coefficient, across the `m` fits in
`long_fit_all_list`:

- the pooled estimate is the mean of the `m` per-completion estimates
  (\\\bar{Q}\\);

- the pooled variance (\\T\\) adds the average within-imputation
  variance (\\\bar{U}\\, from each fit's own standard error) to the
  between-imputation variance (\\B\\, the sample variance of the `m`
  estimates) inflated by a factor of \\1 + 1/m\\: \\T = \bar{U} + (1 +
  1/m) B\\;

- the pooled degrees of freedom use the Barnard & Rubin (1999)
  adjustment, which interpolates between the classical Rubin (1987)
  degrees of freedom (based only on `m` and the fraction of missing
  information) and the complete-data degrees of freedom reported by each
  [`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) fit. This
  adjustment matters most at the small sample sizes typical of this
  package's use case, where the classical, unadjusted degrees of freedom
  can otherwise come out larger than the complete-data degrees of
  freedom, which is not sensible.

A large fraction of missing information (`fmi`, close to 1) for a
coefficient means the missing data cost a lot of precision for that
coefficient specifically, even if the biomarker's overall missingness
rate is modest.

## References

Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys*.
John Wiley & Sons.

Barnard, J. and Rubin, D. B. (1999). Small-Sample Degrees of Freedom
with Multiple Imputation. *Biometrika*, 86(4):948-955.

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
  pooled
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
