# Breslow baseline cumulative hazard for a fixed linear predictor

Breslow baseline cumulative hazard for a fixed linear predictor

## Usage

``` r
breslow_basehaz(time, status, lp)
```

## Arguments

- time, status:

  Follow-up time and 0/1 event indicator.

- lp:

  Linear predictor (uncentered).

## Value

A data frame with columns `hazard` and `time`, at the distinct follow-up
times, as `basehaz(centered = FALSE)` would give.
