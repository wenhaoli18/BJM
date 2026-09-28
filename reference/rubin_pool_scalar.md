# Combine one coefficient's per-completion estimates with Rubin's rules

Internal helper for
[`poolLongitudinalSub`](https://wenhaoli18.github.io/BJM/reference/poolLongitudinalSub.md):
implements the scalar-parameter pooling rules of Rubin (1987, Ch. 3) and
the Barnard & Rubin (1999) small-sample degrees-of-freedom adjustment,
for a single fixed-effect coefficient's `m` per-completion
estimates/standard errors.

## Usage

``` r
rubin_pool_scalar(estimates, std_errors, dfcom)
```

## Arguments

- estimates:

  Length-`m` numeric vector of per-completion point estimates.

- std_errors:

  Length-`m` numeric vector of per-completion standard errors (the
  square root of the within-imputation variance).

- dfcom:

  The complete-data degrees of freedom (assumed common across
  completions; see
  [`poolLongitudinalSub`](https://wenhaoli18.github.io/BJM/reference/poolLongitudinalSub.md)).
