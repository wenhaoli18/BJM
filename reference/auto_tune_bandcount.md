# Auto-select "auto" bandcount arguments by doubling until convergence

Shared implementation backing `bandcount1`/
`bandcount2`/`bandcount3 = "auto"` support in
[`predictRisk()`](https://liwh0904.github.io/BJM/reference/predictRisk.md)/[`dynamicPredictionBio()`](https://liwh0904.github.io/BJM/reference/dynamicPredictionBio.md),
and the bandcount pre-resolution done once, up front, by
[`predictPlot()`](https://liwh0904.github.io/BJM/reference/predictPlot.md)/
[`riskPlot()`](https://liwh0904.github.io/BJM/reference/riskPlot.md) (so
their internal horizon/landmark loops do not repeat the auto-tuning
search on every iteration).

Starts every argument named in `auto_names` at its entry in
[`bandcount_auto_start()`](https://liwh0904.github.io/BJM/reference/bandcount_auto_start.md),
doubles all of them together each round, and compares consecutive
results with
[`max_relative_diff()`](https://liwh0904.github.io/BJM/reference/max_relative_diff.md)
until the largest relative change drops below `tol`, or `max_rounds`
extra doublings have been tried – a hard cap, so this never loops
indefinitely: at most `max_rounds + 1` calls to `predict_fun` (the
default `max_rounds = 2` means at most 3 calls). If the cap is hit
without converging, a warning is issued and the result/bandcount at the
largest value tried is returned anyway, rather than erroring, so
automated pipelines are not interrupted.

## Usage

``` r
auto_tune_bandcount(
  predict_fun,
  args,
  auto_names,
  tol = 0.01,
  max_rounds = 2,
  multiplier = 2
)
```

## Arguments

- predict_fun:

  `predictRisk` or `dynamicPredictionBio`.

- args:

  A named list of all of `predict_fun`'s arguments (typically
  `as.list(environment())` captured right after argument validation,
  before any other local variables are created).

- auto_names:

  Character vector naming which element(s) of `args` to auto-tune (e.g.
  `"bandcount1"`, or `c("bandcount1", "bandcount2")`).

## Value

A list with `result` (`predict_fun`'s return value at the resolved
bandcount) and `bandcount` (a named list of the resolved numeric
bandcount value(s), one per element of `auto_names`).
