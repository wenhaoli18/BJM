# Fit `nlme::lme()`, retrying with more optimizer iterations

`lme()` with `opt = "optim"` stops after `msMaxIter` (default 50)
optimizer iterations and then fails with "optim problem, convergence
error code = 1". That happens when a random-effect variance is close to
zero – e.g. a random slope in `long_sub_random` for data with little
slope variation, more often with few measurements per subject. The call
is evaluated as given first, so a fit that succeeds is exactly what it
was before (call included); only if it fails to converge is it refit
with `msMaxIter = 1000`, with a warning. If that fails too, the error
suggests a simpler random-effects formula. Errors that are not about
convergence are passed on unchanged, without a retry.

## Usage

``` r
lme_with_retry(call_expr, label)
```

## Arguments

- call_expr:

  The `nlme::lme(...)` call (unevaluated; evaluated in the caller's
  frame, as if written there).

- label:

  The biomarker's name, for messages.

## Value

The `lme` fit.
