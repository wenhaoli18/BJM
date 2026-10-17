# Variance-covariance matrix Reference: package "lmm" and package "joineRML" function `mvlme`.

The EM is accelerated by SQUAREM (Varadhan and Roland, 2008), and stops
once the relative change in `D` over an accelerated cycle (three EM
updates) is at most `tol.em`, or after `max.iter` EM updates with a
warning.

## Usage

``` r
longitudinalSubVar(thetaLong, l, tol.em, verbose, max.iter = 1000)
```
