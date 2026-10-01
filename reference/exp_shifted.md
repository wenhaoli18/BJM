# Exponentiate log densities after a per-patient shift

Exponentiate log densities after a per-patient shift

## Usage

``` r
exp_shifted(log_density, shift)
```

## Arguments

- log_density:

  A matrix of log densities, columns = patients.

- shift:

  Per-patient shifts, from
  [`patient_log_shift()`](https://wenhaoli18.github.io/BJM/reference/patient_log_shift.md).

## Value

`exp(log_density - shift)`, column-wise.
