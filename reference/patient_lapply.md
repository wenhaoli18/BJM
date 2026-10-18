# Apply a function to every patient, on one core or several

Used by the copula density functions, whose cost is per patient. With
`getOption("BJM.n_cores")` above 1 (set by the `n_cores` argument of
[`predictRisk()`](https://wenhaoli18.github.io/BJM/reference/predictRisk.md)
and
[`predictLongitudinal()`](https://wenhaoli18.github.io/BJM/reference/predictLongitudinal.md))
the patients are spread over that many cores with
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
(serially on Windows). Each patient is computed with its own fixed
random-number seed, derived from its id
([`patient_seed()`](https://wenhaoli18.github.io/BJM/reference/patient_seed.md)),
and the caller's random-number state is restored afterwards, so the
Monte Carlo box probabilities still used with two or more ordinal random
effects give the same result serially, in parallel, and on every run,
without [`set.seed()`](https://rdrr.io/r/base/Random.html).

## Usage

``` r
patient_lapply(ids, FUN)
```

## Arguments

- ids:

  Patient ids.

- FUN:

  Function of one id.

## Value

A list, one element per id.
