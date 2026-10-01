# Variance-covariance matrix Reference: package "lmm" and package "joineRML" function `mvlme`.

The EM loop stops after `max.iter` iterations even if the relative
change in `D` is still above `tol.em` (previously it had no limit, so a
fit that never met the tolerance never returned), and warns when that
happens.

## Usage

``` r
longitudinalSubVar(thetaLong, l, tol.em, verbose, max.iter = 1000)
```
