# Build a survival-time transform basis from cut points

[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md),
[`dynamicPredictionBio()`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md),
[`predictPlot()`](https://wenhaoli18.github.io/BJM/reference/predictPlot.md),
and
[`riskPlot()`](https://wenhaoli18.github.io/BJM/reference/riskPlot.md)
all take a pair of arguments,
`survival_variable_all`/`survival_trans_function`, that describe
transformed basis variables of the (remaining) survival time; these can
optionally be referenced in `long_sub_fixed` formulas to let the
longitudinal sub-model depend flexibly on time-to-event. In every
example in this package, that pair follows the same convention:
variables named `"Tyears1"`, `"Tyears2"`, ... , each defined as the
absolute distance from a fixed cut point (`function(x) abs(x - k)`).
`survivalTrans()` builds exactly that pair from a plain vector of cut
points, so you do not have to hand-write two parallel lists of matching
names and closures. You remain free to construct
`survival_variable_all`/`survival_trans_function` by hand for any other
transform.

## Usage

``` r
survivalTrans(cut_points, prefix = "Tyears")
```

## Arguments

- cut_points:

  A non-empty numeric vector of cut points, one per transformed basis
  variable.

- prefix:

  Prefix used for the generated variable names in
  `survival_variable_all` (`"Tyears"` by default, giving `"Tyears1"`,
  `"Tyears2"`, ...).

## Value

A named list with elements `survival_variable_all` and
`survival_trans_function`, in the format expected by
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md),
[`dynamicPredictionBio()`](https://wenhaoli18.github.io/BJM/reference/dynamicPredictionBio.md),
[`predictPlot()`](https://wenhaoli18.github.io/BJM/reference/predictPlot.md),
and
[`riskPlot()`](https://wenhaoli18.github.io/BJM/reference/riskPlot.md).

## Examples

``` r
trans <- survivalTrans(c(1, 3, 5, 7))
trans$survival_variable_all
#> [[1]]
#> [1] "Tyears1"
#> 
#> [[2]]
#> [1] "Tyears2"
#> 
#> [[3]]
#> [1] "Tyears3"
#> 
#> [[4]]
#> [1] "Tyears4"
#> 
trans$survival_trans_function[[1]](2)
#> [1] 1

# equivalent to hand-writing:
survival_variable_all <- list("Tyears1", "Tyears2", "Tyears3", "Tyears4")
survival_trans_function <- list(
  fun1 = function(x) abs(x - 1),
  fun2 = function(x) abs(x - 3),
  fun3 = function(x) abs(x - 5),
  fun4 = function(x) abs(x - 7)
)

# survivalTrans() only ever builds the abs(x - k) family above. Any other
# transform -- as long as it is a function of a single numeric time value
# that returns a single finite numeric value -- is written by hand the
# same way, one entry of survival_trans_function per entry of
# survival_variable_all. A few examples that all work equally well with
# predictRisk()/dynamicPredictionBio()/predictPlot()/riskPlot():
survival_variable_all <- list("Tlog", "Tsqrt", "Tsq", "Texp", "Tinv")
survival_trans_function <- list(
  fun1 = function(x) log(x + 1),        # log(x + 1)
  fun2 = function(x) sqrt(abs(x)),      # sqrt(abs(x))
  fun3 = function(x) x^2,               # x^2
  fun4 = function(x) exp(-x / 10),      # exp(-x / 10)
  fun5 = function(x) 1 / (x + 1)        # 1 / (x + 1)
)
sapply(survival_trans_function, function(f) f(2))
#>      fun1      fun2      fun3      fun4      fun5 
#> 1.0986123 1.4142136 4.0000000 0.8187308 0.3333333 
```
