# Summary method for `survivalSub.BJM` objects

Called via `summary(survival_fit_all)`. Returns (and prints) an extended
summary including baseline hazard range, BIC, and McFadden R2 for the
competing-risks GLM.

## Usage

``` r
# S3 method for class 'survivalSub.BJM'
summary(object, digits = 4, ...)
```

## Arguments

- object:

  A `survivalSub.BJM` object returned by
  [`survivalSub`](https://wenhaoli18.github.io/BJM/reference/survivalSub.md).

- digits:

  Number of significant digits. Default is 4.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns a list with components `cox_summary` and (if competing
risks) `glm_summary`.

## Examples

``` r
data(pbc3)
data_survival_fitting <- pbc3[!duplicated(pbc3$id), ]
form_marginal_surv  <- Surv(years, status3) ~ age + sex
form_conditional_cr <- status4 ~ years + age + sex
survival_fit_all  <- survivalSub(data_survival_fitting,
                                 form_marginal_surv, form_conditional_cr)
summary(survival_fit_all)
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
#> <environment: 0x55ab84daa270>
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
#>   Baseline cumulative hazard range: [0.0017, 0.6669]
#>   Time range: [0.1123, 14.3057]
#> 
#> =================================================================
#>  Conditional Competing-Risks Sub-model  [Logistic GLM]
#> -----------------------------------------------------------------
#>  Formula: status4 ~ years + age + sex
#> <environment: 0x55ab84daa270>
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
#>   BIC : 135.57
#>   McFadden R2 : 0.2574
#> =================================================================
#> 
```
