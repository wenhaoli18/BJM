# Summary method for `longitudinalSub.BJM` objects

Like `print` but adds per-outcome random-effects variance components and
the full correlation matrix of D.

## Usage

``` r
# S3 method for class 'longitudinalSub.BJM'
summary(object, digits = 4, ...)
```

## Arguments

- object:

  A `longitudinalSub.BJM` object.

- digits:

  Number of significant digits. Default is 4.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns a list of per-outcome `summary.lme` objects.
