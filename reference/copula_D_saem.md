# Copula EM by stochastic approximation (SAEM)

Fallback of
[`longitudinalSubVarCopula()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSubVarCopula.md)
for many ordinal random effects. Each iteration draws every ordinal
latent score from its truncated normal given the subject's current
random effects (one dimension at a time, as they are independent given
\\b\\), computes \\E\[b b^T \mid y, z\]\\ exactly, draws new random
effects from their Gaussian posterior, and updates \\D\\ by stochastic
approximation: step 1 for `n_burnin` iterations, then \\1/k\\, averaging
the next `n_average`. The draws use a fixed seed, with the caller's
random number state restored, so a fit is reproducible.

## Usage

``` r
copula_D_saem(D, setup, n_burnin = 100, n_average = 200, verbose = FALSE)
```

## Arguments

- D:

  Starting covariance.

- setup:

  Output of
  [`copula_estep_setup()`](https://wenhaoli18.github.io/BJM/reference/copula_estep_setup.md).

- n_burnin, n_average:

  Iterations with step 1 and with step \\1/k\\.

- verbose:

  Print the change in `D` at every iteration.

## Value

A list with `D`.
