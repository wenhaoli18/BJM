# Stratum of each patient for a stratified interval-censored fit

Stratum of each patient for a stratified interval-censored fit

## Usage

``` r
icph_strata(fit, newdata)
```

## Arguments

- fit:

  An `icph.BJM` object.

- newdata:

  One row per patient.

## Value

Character labels in the fit's
[`strata()`](https://rdrr.io/pkg/survival/man/strata.html) format, or
`NULL` when the fit is not stratified.
