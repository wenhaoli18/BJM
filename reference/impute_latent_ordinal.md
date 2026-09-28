# Impute the latent Gaussian score for an ordinal biomarker (inner E-step)

The Gaussian-copula extension represents each ordinal biomarker
observation as a censored/interval-observed draw from an underlying
continuous latent Gaussian variable (the classical Albert & Chib probit
data-augmentation representation): the category actually observed only
tells us the latent score fell in the interval implied by the fitted
thresholds, not its exact value. Before each outer E-step/M-step update
of the shared random-effects covariance matrix `D` (see
[`longitudinalSubVarCopula()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSubVarCopula.md)),
this "inner" E-step plugs in, for every ordinal observation, the mean of
that latent score's truncated-normal conditional distribution given its
observed category and the current fixed-effect and (from the previous
outer E-step) predicted random-effect contribution.

This is a deterministic moment-matching plug-in (not a full
Bayesian/MCEM draw), which is a documented v1 simplification: it is the
mean-field approximation classically used to bootstrap probit/ordinal EM
algorithms, but does not propagate the extra uncertainty a stochastic
draw would.

## Usage

``` r
impute_latent_ordinal(m, l, beta, alpha_m, Eb_prev, beta_idx_m, r_idx_m)
```

## Arguments

- m:

  Index of the ordinal biomarker being imputed.

- l:

  The per-subject design-matrix bookkeeping list built by
  [`longitudinalSubCopula()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSubCopula.md)
  (same shape as the `l` passed to
  [`longitudinalSubVar()`](https://wenhaoli18.github.io/BJM/reference/longitudinalSubVar.md));
  `l$yik[[m]]` holds each subject's observed category codes (`1..K`) for
  biomarker `m`.

- beta:

  Current stacked fixed-effect vector (all biomarkers).

- alpha_m:

  Numeric vector of `K - 1` cumulative-link thresholds for biomarker `m`
  (`lfit[[m]]$alpha`).

- Eb_prev:

  Named list of each subject's posterior mean random-effects vector from
  the previous outer E-step, or `NULL` on the first iteration
  (random-effect contribution is then treated as 0).

- beta_idx_m:

  Integer indices of `beta` belonging to biomarker `m`.

- r_idx_m:

  Integer indices of the stacked random-effects vector belonging to
  biomarker `m`.

## Value

A named list (one element per subject id, matching `l$yik[[m]]`) of
imputed latent-score vectors.
