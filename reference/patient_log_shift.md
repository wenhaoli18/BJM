# Per-patient shift for exponentiating log densities

The conditional-density helpers
([`conditionalYT()`](https://wenhaoli18.github.io/BJM/reference/conditionalYT.md)
and relatives) return log densities, one column per patient. Every
quantity a prediction needs is a ratio of sums of those densities within
one patient, so any per-patient constant cancels: this returns, for each
patient, the largest finite log density across all the matrices given,
to subtract before exponentiating (see
[`exp_shifted()`](https://wenhaoli18.github.io/BJM/reference/exp_shifted.md)).
Without it the densities – whose scale is the determinant of a
covariance matrix that grows with the number of observations and with
the biomarkers' units – underflowed to 0 or overflowed to `Inf`.

## Usage

``` r
patient_log_shift(...)
```

## Arguments

- ...:

  Matrices of log densities, rows = grid points, columns = patients (all
  with the same columns).

## Value

A numeric vector, one shift per patient (`0` for a patient with no
finite value).
