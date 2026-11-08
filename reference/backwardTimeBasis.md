# Add a backward-time basis to longitudinal sub-model formulas

The backward joint model describes a biomarker given the event time
\\T\\.
[`survivalTrans()`](https://wenhaoli18.github.io/BJM/reference/survivalTrans.md)
adds flexibility in \\T\\ itself, through terms \\\|T - c\|\\ that are
the same at every visit of a patient. Biomarkers often change instead
with the time left until the event, the *backward time* \\u = T - s\\,
where \\s\\ is the visit time: for example, a marker that stays flat for
years and then accelerates in the last year or two before the event. A
term in \\T\\, \\s\\ and \\T \times s\\ can only describe such a
trajectory as linear in \\u\\.

`backwardTimeBasis()` adds the piecewise-linear basis \\\|u - c\| =
\|T - s - c\|\\, one term per cut point \\c\\, to each formula. The
terms are written as `I(abs(years - year - c))`, so they are ordinary
columns of the fixed-effects design. They need no
`survival_variable_all`/`survival_trans_function` entries:
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md),
[`predictLongitudinal()`](https://wenhaoli18.github.io/BJM/reference/predictLongitudinal.md)
and
[`simulateTrajectory()`](https://wenhaoli18.github.io/BJM/reference/simulateTrajectory.md)
evaluate them at every integration grid point of the event time from the
formula itself, using each row's visit time.

## Usage

``` r
backwardTimeBasis(
  long_sub_fixed,
  cut_points,
  survival_variable = "years",
  time_variable = "year"
)
```

## Arguments

- long_sub_fixed:

  A fixed-effects formula, or a list of them (one per biomarker, as
  passed to
  [`longitudinalSub()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSub.md)).
  Each formula should already use the event time, e.g.
  `serBilir ~ year + age + sex + years`.

- cut_points:

  A non-empty numeric vector of cut points on the backward-time scale,
  in the units of `survival_variable`. Knots where the trajectory may
  bend, e.g. `c(0.5, 1, 2, 4)` for changes concentrated in the last few
  years before the event.

- survival_variable:

  Name of the event-time column (`"years"` in `pbc3`).

- time_variable:

  Name of the visit-time column (`"year"` in `pbc3`).

## Value

`long_sub_fixed` with the terms
`I(abs(survival_variable - time_variable - c))` appended to the
right-hand side of each formula. A single formula in gives a single
formula out, and a list keeps its names. Each formula keeps its
environment.

## See also

[`survivalTrans()`](https://wenhaoli18.github.io/BJM/reference/survivalTrans.md)
for a basis in the event time itself;
[`vignette("nonlinear-terms", package = "BJM")`](https://wenhaoli18.github.io/BJM/articles/nonlinear-terms.md).

## Examples

``` r
backwardTimeBasis(serBilir ~ year + age + sex + years, c(1, 3))
#> serBilir ~ year + age + sex + years + I(abs(years - year - 1)) + 
#>     I(abs(years - year - 3))
#> <environment: 0x62611b6048f0>

long_sub_fixed <- backwardTimeBasis(
  list("serBilir" = serBilir ~ year + age + sex + years,
       "albumin"  = albumin  ~ year + age + sex + years),
  cut_points = c(1, 3, 5)
)
long_sub_fixed
#> $serBilir
#> serBilir ~ year + age + sex + years + I(abs(years - year - 1)) + 
#>     I(abs(years - year - 3)) + I(abs(years - year - 5))
#> <environment: 0x62611b6048f0>
#> 
#> $albumin
#> albumin ~ year + age + sex + years + I(abs(years - year - 1)) + 
#>     I(abs(years - year - 3)) + I(abs(years - year - 5))
#> <environment: 0x62611b6048f0>
#> 

# \donttest{
data(pbc3)
data_fit_all <- list(pbc3[pbc3$status3 == 1, ], pbc3[pbc3$status3 == 1, ])
long_fit_all <- longitudinalSub(data_fit_all, long_sub_fixed,
                                list("serBilir" = ~ year | id, "albumin" = ~ year | id))
# }
```
