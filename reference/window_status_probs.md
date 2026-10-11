# Probability of a case and of a control from the posterior of T

Model-based case/control status of one subject in the window \\(s,
e\]\\, given what was observed: the event in \\(L, R\]\\ with type
`cause`, or no event by \\L\\ (`R = Inf`). With \\F_k\\ the posterior
cumulative mass of event type \\k\\
([`event_time_posterior()`](https://wenhaoli18.github.io/BJM/reference/event_time_posterior.md),
linear within cells), the probability of being a case of type `k` is the
posterior mass of \\(\max(L, s), \min(R, e)\]\\ for the observed type
over that of \\(L, R\]\\, and of being a control (event-free at \\e\\)
the mass beyond \\e\\ over that of \\(L, R\]\\. Yang et al. (2026),
Biometrical Journal 68:e70108.

## Usage

``` r
window_status_probs(post, L, R, cause, s, e, k)
```

## Arguments

- post:

  One element of
  [`event_time_posterior()`](https://wenhaoli18.github.io/BJM/reference/event_time_posterior.md).

- L, R, cause:

  The subject's observed outcome.

- s, e:

  The window.

- k:

  The event type evaluated.

## Value

`c(case = , control = )`.
