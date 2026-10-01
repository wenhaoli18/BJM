# Check whether bandcount1/bandcount2/bandcount3 are large enough

[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)/[`predictLongitudinal()`](https://wenhaoli18.github.io/BJM/reference/predictLongitudinal.md)
approximate the integrals behind the predicted risk probabilities (and,
for
[`predictLongitudinal()`](https://wenhaoli18.github.io/BJM/reference/predictLongitudinal.md),
the predicted biomarker density) with a finite grid, controlled by
`bandcount1`/`bandcount2`/ `bandcount3`. There is no universal correct
value: too few grid points silently bias the answer, and too many just
cost more time, and the right value depends on the data (e.g. how wide
the follow-up range is). Rather than guess, or auto-loop until some
tolerance is met – which multiplies runtime unpredictably, especially
for
[`predictLongitudinal()`](https://wenhaoli18.github.io/BJM/reference/predictLongitudinal.md)'s
nested `bandcount2` x `bandcount3` grid – this runs the prediction once
at the bandcount value(s) you supply, once more with those value(s)
scaled up, and reports the largest relative change between the two, so
you can see whether you have already converged or need to increase the
checked bandcount(s) and re-run. It costs exactly 2 calls to
`predict_fun`, regardless of how many times you invoke it.

## Usage

``` r
checkBandcountConvergence(
  predict_fun,
  ...,
  bandcount_args,
  multiplier = 2,
  tol = 0.01
)
```

## Arguments

- predict_fun:

  The prediction function to check: `predictRisk` or
  `predictLongitudinal` themselves (not a string, and not
  [`predictPlot()`](https://wenhaoli18.github.io/BJM/reference/predictPlot.md)/[`riskPlot()`](https://wenhaoli18.github.io/BJM/reference/riskPlot.md),
  which return a plot rather than the underlying numeric predictions
  this function compares). When passing `predictLongitudinal`, forward a
  `bio_i` (via `...`) that names exactly **one** biomarker – this
  function only supports comparing a single biomarker's prediction at a
  time (see `...`).

- ...:

  Arguments to forward to `predict_fun`, exactly as you would call it
  directly, except for the bandcount argument(s) being checked, which
  are supplied separately via `bandcount_args`. Any bandcount argument
  of `predict_fun` not named in `bandcount_args` must still be supplied
  here (it is held fixed at that value for both calls).

- bandcount_args:

  A named list giving the bandcount value(s) to check, e.g.
  `list(bandcount1 = 10, bandcount2 = 40)`. Every element is scaled by
  `multiplier` for the second call. To isolate which bandcount is
  driving instability, check one at a time (e.g.
  `list(bandcount2 = 40)`, with `bandcount3` passed as a fixed value via
  `...`), rather than checking all of them together.

- multiplier:

  Factor the checked bandcount value(s) are scaled by for the second
  call. Must be greater than 1. Defaults to `2` (doubling).

- tol:

  Relative-change tolerance below which the result is reported as
  converged. Defaults to `0.01` (1%).

## Value

An object of class `"checkBandcountConvergence.BJM"`, a named list with
elements:

- base_bandcount:

  The bandcount value(s) checked, as supplied in `bandcount_args`.

- scaled_bandcount:

  `base_bandcount` scaled by `multiplier`.

- tol:

  The relative-change tolerance used to decide `converged`.

- by_field:

  A named numeric vector giving the largest relative change, per
  comparable output field (e.g. `risk_prob_1`, `Y_predict`), between the
  base and scaled call.

- max_rel_diff:

  The largest value in `by_field`, or `NA` if there was no comparable
  output (e.g. `horizon <= 0`).

- converged:

  `TRUE` if `max_rel_diff < tol`, `FALSE` if not, `NA` if there was
  nothing to compare.

- base_result:

  The full result of calling `predict_fun` at `base_bandcount`.

- scaled_result:

  The full result of calling `predict_fun` at `scaled_bandcount`.

Printing the object summarizes these fields.

## Examples

``` r

# \donttest{
data(pbc3)

data_survival_fitting = pbc3[!duplicated(pbc3$id), ]
survival_fit_all = survivalSub(data_survival_fitting, Surv(years, status3) ~ age + sex, NULL)

long_sub_fixed = serBilir ~ year + age + sex + (years) + (years) * year
long_sub_random = ~ year | id
long_fit_all = longitudinalSub(pbc3[pbc3$status3 == 1, ], long_sub_fixed, long_sub_random)

trans = survivalTrans(c(1, 3, 5, 7))

data_predict_all = pbc3[pbc3$id == 2 & pbc3$year <= 3, ]

# Check bandcount1/bandcount2 together: are they both already large enough?
check = checkBandcountConvergence(
  predictRisk, data_predict_all, long_fit_all, survival_fit_all,
  prediction_time = 3, horizon = 3, time_variable = "year",
  trans$survival_variable_all, trans$survival_trans_function,
  bandcount_args = list(bandcount1 = 10, bandcount2 = 10)
)
print(check)
#> Bandcount convergence check
#>   base:   bandcount1 = 10, bandcount2 = 10
#>   scaled: bandcount1 = 20, bandcount2 = 20
#>   max relative change in risk_prob_1: 0.0472
#>   NOT converged: max relative change (0.0472) exceeds tol (0.01); consider increasing the checked bandcount(s) and re-running.
# }
```
