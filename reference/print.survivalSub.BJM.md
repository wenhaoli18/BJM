# Print method for `survivalSub.BJM` objects

Automatically called when you type `survival_fit_all` or
`print(survival_fit_all)` at the console. Displays a JMbayes2-style
formatted summary of the survival sub-model.

## Usage

``` r
# S3 method for class 'survivalSub.BJM'
print(x, digits = 4, ...)
```

## Arguments

- x:

  A `survivalSub.BJM` object returned by
  [`survivalSub`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md).

- digits:

  Number of significant digits. Default is 4.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns `x`.

## Examples

``` r
data(pbc3)
data_survival_fitting <- pbc3[!duplicated(pbc3$id), ]
form_marginal_surv  <- Surv(years, status3) ~ age + sex
form_conditional_cr <- status4 ~ years + age + sex
survival_fit_all  <- survivalSub(data_survival_fitting,
                                 form_marginal_surv, form_conditional_cr)
survival_fit_all   # triggers print.survivalSub.BJM automatically
#> 
#> Call:
#> survivalSub(form_marginal_surv = Surv(years, status3) ~ age + sex,
#>             form_conditional_cr = status4 ~ years + age + sex)
#> 
#> Data Descriptives:
#>   Number of subjects        : 312
#>   Number of events          : 169
#>   Cause-1 events (CR model) : 29
#>   Cause-2 events (CR model) : 140
#> 
#> =================================================================
#>  Marginal Survival Sub-model  [Cox PH]
#> -----------------------------------------------------------------
#>  Formula: Surv(years, status3) ~ age + sex
#> <environment: 0x5620c5494f00>
#> 
#>          Coef exp(Coef)        SE      z p-value   
#> age  0.020411  1.020621  0.007584  2.691 0.00712 **
#> sex -0.497707  0.607923  0.207078 -2.403 0.01624 * 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>   n = 312,  events = 169
#>   Concordance       = 0.570  (se = 0.0226)
#>   Likelihood ratio  = 14.20  on 2 df,  p = 0.0008263
#>   Wald test         = 15.29  on 2 df,  p = 0.0004776
#>   Score (logrank)   = 15.50  on 2 df,  p = 0.0004299
#> 
#> =================================================================
#>  Conditional Competing-Risks Sub-model  [Logistic GLM]
#> -----------------------------------------------------------------
#>  Formula: status4 ~ years + age + sex
#> <environment: 0x5620c5494f00>
#> 
#>                 Coef       SE      z  p-value    
#> (Intercept)  5.65622  1.70478  3.318 0.000907 ***
#> years       -0.02898  0.08919 -0.325 0.745242    
#> age         -0.15305  0.03076 -4.976  6.5e-07 ***
#> sex          0.02109  0.77150  0.027 0.978192    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>   Null deviance     : 154.94  on 168 df
#>   Residual deviance : 115.05  on 165 df
#>   AIC : 123.05
#> =================================================================
#> 
```
